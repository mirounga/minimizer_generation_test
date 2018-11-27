#include "stdafx.h"
#include "CppUnitTest.h"

#include"stdio.h"

#include <algorithm>
#include <chrono>
#include <vector>
#include <deque>

#include "intrin.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace minimizer_test
{
	void ReadSequence(std::vector<char> &seq, const char* fname)
	{
		FILE* fseq;
		fopen_s(&fseq, fname, "rb");

		fseek(fseq, 0, SEEK_END);
		std::size_t seq_len = ftell(fseq);
		fseek(fseq, 0, SEEK_SET);

		seq.resize(seq_len);
		fread_s(seq.data(), seq_len, sizeof(char), seq_len, fseq);
	}

	inline uint32_t get_seed(const char* seq, uint64_t kmer_mask, std::size_t pos)
	{
		uint64_t seed_qword = *((uint64_t*)(seq + (pos >> 2)));
		uint32_t seed = (uint32_t)((seed_qword >> ((pos & 0x03ull) << 1)) & kmer_mask);

		// hash32
		seed = (~seed + (seed << 21)) & kmer_mask;
		seed = seed ^ (seed >> 24);
		seed = ((seed + (seed << 3)) + (seed << 8)) & kmer_mask;
		seed = seed ^ (seed >> 14);
		seed = ((seed + (seed << 2)) + (seed << 4)) & kmer_mask;
		seed = seed ^ (seed >> 28);
		seed = (seed + (seed << 31)) & kmer_mask;

		return seed;
	}

	std::size_t generate_minimizers(std::vector<char> &seq, std::size_t kmer_size, std::size_t window_size, uint32_t * minimizers)
	{
		uint32_t last_m = 0;
		uint64_t last_p = 0;

		uint64_t kmer_mask = (1ull << 2ull * kmer_size) - 1ull;

		uint32_t* window = (uint32_t*)calloc(window_size, sizeof(uint32_t));

		std::size_t num_min = 0, p = 0;

		for (p = 0; p < window_size - 1; p++)
		{
			window[p] = get_seed(seq.data(), kmer_mask, p);
		}

		for (; p < 4 * seq.size() - kmer_size; p++)
		{
			window[p%window_size] = get_seed(seq.data(), kmer_mask, p);

			uint32_t m = ~0ul;
			for (int i = 0; i < window_size; i++)
			{
				if (m > window[i])
				{
					m = window[i];
				}
			}

			if ((m != last_m) || (p - last_p >= window_size))
			{
				minimizers[num_min++] = m;
				last_m = m;
				last_p = p;
			}
		}

		free(window);

		return num_min;
	}

	std::size_t generate_minimizers_q(std::vector<char> &seq, std::size_t kmer_size, std::size_t window_size, uint32_t * minimizers)
	{
		uint32_t last_m = 0;
		uint64_t last_p = 0;

		uint64_t kmer_mask = (1ull << 2ull * kmer_size) - 1ull;

		std::deque<std::pair<uint32_t, uint64_t>> min_window;

		std::size_t num_min = 0;

		for (std::size_t p = 0; p < 4 * seq.size() - kmer_size; p++)
		{
			uint32_t seed = get_seed(seq.data(), kmer_mask, p);

			while (!min_window.empty() && min_window.back().first >= seed)
			{
				min_window.pop_back();
			}

			min_window.push_back(std::make_pair(seed, p));

			while (p - min_window.front().second >= window_size)
			{
				min_window.pop_front();
			}

			uint32_t m = min_window.front().first;

			if ((p >= window_size - 1) && ((m != last_m) || (p - last_p >= window_size)))
			{
				minimizers[num_min++] = m;
				last_m = m;
				last_p = p;
			}
		}

		return num_min;
	}

#define _mm256_shuffle_epi32(a, b, imm) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), imm))

	__forceinline __m256i _get_seeds(const __m256i _seed_qword, const __m256i _kmer_mask)
	{
		// Vector of 8 seeds
		__m256i _seed_x_8 = _mm256_and_si256(
			_mm256_shuffle_epi32(
				_mm256_srlv_epi64(
					_seed_qword,
					_mm256_set_epi64x(
						10ull, 8ull, 2ull, 0ull)),
				_mm256_srlv_epi64(
					_seed_qword,
					_mm256_set_epi64x(
						14ull, 12ull, 6ull, 4ull)),
				0x88),
			_kmer_mask);

		// hash32
		//seed = (~seed + (seed << 21)) & m;
		_seed_x_8 = _mm256_and_si256(
			_mm256_add_epi32(
				_mm256_xor_si256(
					_seed_x_8,
					_mm256_set1_epi32(-1)),
				_mm256_slli_epi32(
					_seed_x_8,
					21)),
			_kmer_mask);

		//seed = seed ^ (seed >> 24);
		_seed_x_8 = _mm256_xor_si256(
			_seed_x_8,
			_mm256_srli_epi32(
				_seed_x_8,
				24));

		//seed = ((seed + (seed << 3)) + (seed << 8)) & m;
		_seed_x_8 = _mm256_and_si256(
			_mm256_add_epi32(
				_seed_x_8,
				_mm256_add_epi32(
					_mm256_slli_epi32(
						_seed_x_8,
						3),
					_mm256_slli_epi32(
						_seed_x_8,
						8))),
			_kmer_mask);

		//seed = seed ^ (seed >> 14);
		_seed_x_8 = _mm256_xor_si256(
			_seed_x_8,
			_mm256_srli_epi32(
				_seed_x_8,
				14));

		//seed = ((seed + (seed << 2)) + (seed << 4)) & m;
		_seed_x_8 = _mm256_and_si256(
			_mm256_add_epi32(
				_seed_x_8,
				_mm256_add_epi32(
					_mm256_slli_epi32(
						_seed_x_8,
						2),
					_mm256_slli_epi32(
						_seed_x_8,
						4))),
			_kmer_mask);

		//seed = seed ^ (seed >> 28);
		_seed_x_8 = _mm256_xor_si256(
			_seed_x_8,
			_mm256_srli_epi32(
				_seed_x_8,
				28));

		//seed = (seed + (seed << 31)) & m;
		_seed_x_8 = _mm256_and_si256(
			_mm256_add_epi32(
				_seed_x_8,
				_mm256_slli_epi32(
					_seed_x_8,
					31)),
			_kmer_mask);

		return _seed_x_8;
	}

	__forceinline __m256i _sliding_minimum_5w(__m256i &_min_window, const __m256i _min_seeds)
	{
		static const __m256i _mm256_000f_epu32 = _mm256_set_epi32(
			0, 0, 0, 0xffffffff,
			0, 0, 0, 0xffffffff);

		// Suffix scan
		__m256i _min_suffix = _mm256_min_epu32(
			_min_seeds,
			_mm256_or_si256(
				_mm256_bsrli_epi128(
					_min_seeds,
					4),
				_mm256_bslli_epi128(
					_mm256_000f_epu32,
					12)));

		_min_suffix = _mm256_min_epu32(
			_min_suffix,
			_mm256_shuffle_epi32(
				_min_suffix,
				_mm256_000f_epu32,
				0x0e));

		// Prefix scan
		__m256i _min_prefix = _mm256_min_epu32(
			_mm256_or_si256(
				_mm256_bslli_epi128(
					_min_seeds,
					4),
				_mm256_000f_epu32),
			_min_seeds);

		_min_prefix = _mm256_min_epu32(
			_mm256_shuffle_epi32(
				_mm256_000f_epu32,
				_min_prefix,
				0x40),
			_min_prefix);

		// Sliding window minimum
		__m256i _min_result = _mm256_min_epu32(
			_mm256_permute2x128_si256(
				_min_window,
				_min_suffix,
				0x21),
			_min_prefix);

		_min_window = _min_suffix;

		return _min_result;
	}

	template <size_t __N>
	__forceinline __m256i _mm256_shift_left_si256(__m256i a, __m256i b) {
		__m256i c = _mm256_permute2x128_si256(a, b, 0x03);
		return _mm256_alignr_epi8(a, c, 16 - __N);
	}

	template <size_t __N>
	__forceinline __m256i _mm256_shift_right_si256(__m256i a, __m256i b) {
		__m256i c = _mm256_permute2x128_si256(a, b, 0x21);
		return _mm256_alignr_epi8(c, a, __N);
	}

	inline __m256i _sliding_minimum_9w(__m256i &_min_window, const __m256i _min_seeds)
	{
		static const __m256i _mm256_ff_epu32 = _mm256_set1_epi32(0xffffffff);

		// Suffix scan
		__m256i _min_suffix = _mm256_min_epu32(
			_min_seeds,
			_mm256_shift_right_si256<4>(
				_min_seeds,
				_mm256_ff_epu32));

		_min_suffix = _mm256_min_epu32(
			_min_suffix,
			_mm256_shift_right_si256<8>(
				_min_suffix,
				_mm256_ff_epu32));

		_min_suffix = _mm256_min_epu32(
			_min_suffix,
			_mm256_permute2x128_si256(
				_min_suffix,
				_mm256_ff_epu32,
				0x21));

		// Prefix scan
		__m256i _min_prefix = _mm256_min_epu32(
			_min_seeds,
			_mm256_shift_left_si256<4>(
				_min_seeds,
				_mm256_ff_epu32));

		_min_prefix = _mm256_min_epu32(
			_min_prefix,
			_mm256_shift_left_si256<8>(
				_min_prefix,
				_mm256_ff_epu32));

		_min_prefix = _mm256_min_epu32(
			_min_prefix,
			_mm256_permute2x128_si256(
				_min_prefix,
				_mm256_ff_epu32,
				0x03));

		// Sliding window minimum
		__m256i _min_result = _mm256_min_epu32(
			_min_window,
			_min_prefix);

		_min_window = _min_suffix;

		return _min_result;
	}

	inline void _sliding_minimum_17w(__m256i &_min_window0, __m256i &_min_window1, const __m256i _min_seeds0, const __m256i _min_seeds1, __m256i &_min_result0, __m256i &_min_result1)
	{
		static const __m256i _mm256_ff_epu32 = _mm256_set1_epi32(0xffffffff);

		// Suffix scan
		__m256i _min_suffix0 = _mm256_min_epu32(
			_min_seeds0,
			_mm256_shift_right_si256<4>(
				_min_seeds0,
				_min_seeds1));

		__m256i _min_suffix1 = _mm256_min_epu32(
			_min_seeds1,
			_mm256_shift_right_si256<4>(
				_min_seeds1,
				_mm256_ff_epu32));

		_min_suffix0 = _mm256_min_epu32(
			_min_suffix0,
			_mm256_shift_right_si256<8>(
				_min_suffix0,
				_min_suffix1));

		_min_suffix1 = _mm256_min_epu32(
			_min_suffix1,
			_mm256_shift_right_si256<8>(
				_min_suffix1,
				_mm256_ff_epu32));

		_min_suffix0 = _mm256_min_epu32(
			_min_suffix0,
			_mm256_permute2x128_si256(
				_min_suffix0,
				_min_suffix1,
				0x21));

		_min_suffix1 = _mm256_min_epu32(
			_min_suffix1,
			_mm256_permute2x128_si256(
				_min_suffix1,
				_mm256_ff_epu32,
				0x21));

		_min_suffix0 = _mm256_min_epu32(
			_min_suffix0,
			_min_suffix1);

		// Prefix scan
		__m256i _min_prefix1 = _mm256_min_epu32(
			_min_seeds1,
			_mm256_shift_left_si256<4>(
				_min_seeds1,
				_min_seeds0));

		__m256i _min_prefix0 = _mm256_min_epu32(
			_min_seeds0,
			_mm256_shift_left_si256<4>(
				_min_seeds0,
				_mm256_ff_epu32));

		_min_prefix1 = _mm256_min_epu32(
			_min_prefix1,
			_mm256_shift_left_si256<8>(
				_min_prefix1,
				_min_prefix0));

		_min_prefix0 = _mm256_min_epu32(
			_min_prefix0,
			_mm256_shift_left_si256<8>(
				_min_prefix0,
				_mm256_ff_epu32));

		_min_prefix1 = _mm256_min_epu32(
			_min_prefix1,
			_mm256_permute2x128_si256(
				_min_prefix1,
				_min_prefix0,
				0x03));

		_min_prefix0 = _mm256_min_epu32(
			_min_prefix0,
			_mm256_permute2x128_si256(
				_min_prefix0,
				_mm256_ff_epu32,
				0x03));

		_min_prefix1 = _mm256_min_epu32(
			_min_prefix1,
			_min_prefix0);

		// Sliding window minimum
		_min_result0 = _mm256_min_epu32(
			_min_window0,
			_min_prefix0);

		_min_result1 = _mm256_min_epu32(
			_min_window1,
			_min_prefix1);

		_min_window0 = _min_suffix0;
		_min_window1 = _min_suffix1;
	}

	std::size_t generate_minimizers_qw(std::vector<char> &seq, std::size_t kmer_size, std::size_t window_size, uint32_t * minimizers)
	{
		uint32_t last_m = 0;
		uint64_t last_p = 0;

		std::deque<std::pair<uint32_t, uint64_t>> min_window;

		std::size_t num_min = 0;

		uint64_t qlen_centinel = (4 * seq.size()) - kmer_size;

		__m256i _kmer_mask = _mm256_set1_epi32(
			(1 << (kmer_size << 1)) - 1);

		// Main loop
		for (uint64_t p = 0; p < qlen_centinel;)
		{
			uint64_t seed_qword = *((uint64_t*)(seq.data() + (p >> 2)));

			__m256i _seed_qword = _mm256_set1_epi64x(seed_qword);

			__m256i _min_seeds = _get_seeds(_seed_qword, _kmer_mask);

			for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
			{
				uint32_t seed = _min_seeds.m256i_u32[i];

				while (!min_window.empty() && min_window.back().first >= seed)
				{
					min_window.pop_back();
				}

				min_window.push_back(std::make_pair(seed, p));

				while (p - min_window.front().second >= window_size)
				{
					min_window.pop_front();
				}

				uint32_t m = min_window.front().first;

				if ((p >= window_size - 1) && ((m != last_m) || (p - last_p >= window_size)))
				{
					minimizers[num_min++] = m;
					last_m = m;
					last_p = p;
				}
			}
		}

		return num_min;
	}

	// Special case with w set to 5
	std::size_t generate_minimizers_5w(std::vector<char> &seq, uint64_t kmer_size, uint32_t * minimizers)
	{
		const int64_t window_size = 5ull;

		uint32_t last_m = 0;
		uint64_t last_p = 0;

		std::size_t num_min = 0;

		uint64_t qlen_centinel = (4 * seq.size()) - kmer_size;

		__m256i _kmer_mask = _mm256_set1_epi32(
			(1 << (kmer_size << 1)) - 1);

		__m256i _min_window = _mm256_setzero_si256();

		// Main loop
		for (uint64_t p = 0; p < qlen_centinel;)
		{
			uint64_t seed_qword = *((uint64_t*)(seq.data() + (p >> 2)));

			__m256i _seed_qword = _mm256_set1_epi64x(seed_qword);

			__m256i _min_seeds = _get_seeds(_seed_qword, _kmer_mask);

			__m256i _min_result = _sliding_minimum_5w(_min_window, _min_seeds);

			for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
			{
				uint32_t m = _min_result.m256i_u32[i];

				if ((m != last_m) || (p - last_p >= window_size))
				{
					minimizers[num_min++] = m;
					last_m = m;
					last_p = p;
				}
			}
		}

		return num_min;
	}

	// Special case with w set to 9
	std::size_t generate_minimizers_9w(std::vector<char> &seq, uint64_t kmer_size, uint32_t * minimizers)
	{
		const int64_t window_size = 9ull;

		uint32_t last_m = 0;
		uint64_t last_p = 0;

		std::size_t num_min = 0;

		uint64_t qlen_centinel = (4 * seq.size()) - kmer_size;

		__m256i _kmer_mask = _mm256_set1_epi32(
			(1 << (kmer_size << 1)) - 1);

		__m256i _min_window = _mm256_setzero_si256();

		// Main loop
		for (uint64_t p = 0; p < qlen_centinel;)
		{
			uint64_t seed_qword = *((uint64_t*)(seq.data() + (p >> 2)));

			__m256i _seed_qword = _mm256_set1_epi64x(seed_qword);

			__m256i _min_seeds = _get_seeds(_seed_qword, _kmer_mask);

			__m256i _min_result = _sliding_minimum_9w(_min_window, _min_seeds);

			for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
			{
				uint32_t m = _min_result.m256i_u32[i];

				if ((m != last_m) || (p - last_p >= window_size))
				{
					minimizers[num_min++] = m;
					last_m = m;
					last_p = p;
				}
			}
		}

		return num_min;
	}

	// Special case with w set to 17
	std::size_t generate_minimizers_17w(std::vector<char> &seq, uint64_t kmer_size, uint32_t * minimizers)
	{
		const int64_t window_size = 17ull;

		uint32_t last_m = 0;
		uint64_t last_p = 0;

		std::size_t num_min = 0;

		uint64_t qlen_centinel = (4 * seq.size()) - kmer_size;

		__m256i _kmer_mask = _mm256_set1_epi32(
			(1 << (kmer_size << 1)) - 1);

		__m256i _min_window0 = _mm256_setzero_si256();
		__m256i _min_window1 = _mm256_setzero_si256();

		// Main loop
		for (uint64_t p = 0; p < qlen_centinel;)
		{
			uint64_t seed_qword = *((uint64_t*)(seq.data() + (p >> 2)));

			__m256i _seed_qword = _mm256_set1_epi64x(seed_qword);

			__m256i _min_seeds0 = _get_seeds(_seed_qword, _kmer_mask);

			seed_qword = *((uint64_t*)(seq.data() + ((p + 8) >> 2)));

			_seed_qword = _mm256_set1_epi64x(seed_qword);

			__m256i _min_seeds1 = _get_seeds(_seed_qword, _kmer_mask);

			__m256i _min_result0, _min_result1;

			_sliding_minimum_17w(_min_window0, _min_window1, _min_seeds0, _min_seeds1, _min_result0, _min_result1);

			for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
			{
				uint32_t m = _min_result0.m256i_u32[i];

				if ((m != last_m) || (p - last_p >= window_size))
				{
					minimizers[num_min++] = m;
					last_m = m;
					last_p = p;
				}
			}

			for (uint32_t i = 0; (i < 8) && (p < qlen_centinel); i++, p++)
			{
				uint32_t m = _min_result1.m256i_u32[i];

				if ((m != last_m) || (p - last_p >= window_size))
				{
					minimizers[num_min++] = m;
					last_m = m;
					last_p = p;
				}
			}
		}

		return num_min;
	}

	TEST_CLASS(VectorizationTest)
	{
		std::size_t test_generic(std::size_t window_size)
		{
			std::vector<char> seq;
			ReadSequence(seq, "data\\b37m1.fa.pac");
//			seq.resize(1000);

			std::vector<uint32_t> minimizers1;
			minimizers1.resize(4 * seq.size());
			auto start = std::chrono::high_resolution_clock::now();

			std::size_t num_min1 = generate_minimizers(seq, 15, window_size, minimizers1.data());

			auto end = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration<double>(end - start).count();

			{
				std::stringstream message;
				message << "Minimizers found: " << num_min1 << std::endl;
				Logger::WriteMessage(message.str().c_str());
			}

			{
				std::stringstream message;
				message << "Array version: " << duration << "s" << std::endl;
				Logger::WriteMessage(message.str().c_str());
			}

			std::vector<uint32_t> minimizers2;
			minimizers2.resize(4 * seq.size());

			start = std::chrono::high_resolution_clock::now();

			std::size_t num_min2 = generate_minimizers_q(seq, 15, window_size, minimizers2.data());

			end = std::chrono::high_resolution_clock::now();

			duration = std::chrono::duration<double>(end - start).count();

			Assert::AreEqual(num_min1, num_min2);

			{
				std::stringstream message;
				message << "Deque version: " << duration << "s" << std::endl;
				Logger::WriteMessage(message.str().c_str());
			}


			std::vector<uint32_t> minimizers3;
			minimizers3.resize(4 * seq.size());

			start = std::chrono::high_resolution_clock::now();

			std::size_t num_min3 = generate_minimizers_qw(seq, 15, window_size, minimizers3.data());

			end = std::chrono::high_resolution_clock::now();

			duration = std::chrono::duration<double>(end - start).count();

			Assert::AreEqual(num_min1, num_min3);

			{
				std::stringstream message;
				message << "Deque + vector version: " << duration << "s" << std::endl;
				Logger::WriteMessage(message.str().c_str());
			}

			return num_min1;
		}

	public:
		TEST_METHOD(TestWindowSize4)
		{
			test_generic(4ull);
		}

		TEST_METHOD(TestWindowSize5)
		{
			std::size_t num_min1 = test_generic(5);

			std::vector<char> seq;
			ReadSequence(seq, "data\\b37m1.fa.pac");

			std::vector<uint32_t> minimizers4;
			minimizers4.resize(num_min1);

			auto start = std::chrono::high_resolution_clock::now();

			std::size_t num_min4 = generate_minimizers_5w(seq, 15, minimizers4.data());

			auto end = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration<double>(end - start).count();

			Assert::AreEqual(num_min1, num_min4);

			{
				std::stringstream message;
				message << "Vector version: " << duration << "s" << std::endl;
				Logger::WriteMessage(message.str().c_str());
			}
		}

		TEST_METHOD(TestWindowSize6)
		{
			test_generic(6ull);
		}

		TEST_METHOD(TestWindowSize7)
		{
			test_generic(7ull);
		}

		TEST_METHOD(TestWindowSize8)
		{
			test_generic(8ull);
		}

		TEST_METHOD(TestWindowSize9)
		{
			std::size_t num_min1 = test_generic(9ull);

			std::vector<char> seq;
			ReadSequence(seq, "data\\b37m1.fa.pac");

			std::vector<uint32_t> minimizers4;
			minimizers4.resize(num_min1);

			auto start = std::chrono::high_resolution_clock::now();

			std::size_t num_min4 = generate_minimizers_9w(seq, 15, minimizers4.data());

			auto end = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration<double>(end - start).count();

			Assert::AreEqual(num_min1, num_min4);

			{
				std::stringstream message;
				message << "Vector version: " << duration << "s" << std::endl;
				Logger::WriteMessage(message.str().c_str());
			}
		}

		TEST_METHOD(TestWindowSize10)
		{
			test_generic(10ull);
		}

		TEST_METHOD(TestWindowSize11)
		{
			test_generic(11ull);
		}

		TEST_METHOD(TestWindowSize12)
		{
			test_generic(12ull);
		}

		TEST_METHOD(TestWindowSize13)
		{
			test_generic(13ull);
		}

		TEST_METHOD(TestWindowSize14)
		{
			test_generic(14ull);
		}

		TEST_METHOD(TestWindowSize15)
		{
			test_generic(15ull);
		}

		TEST_METHOD(TestWindowSize16)
		{
			test_generic(16ull);
		}

		TEST_METHOD(TestWindowSize17)
		{
			std::size_t num_min1 = test_generic(17ull);

			std::vector<char> seq;
			ReadSequence(seq, "data\\b37m1.fa.pac");
//			seq.resize(1000);

			std::vector<uint32_t> minimizers4;
			minimizers4.resize(num_min1);

			auto start = std::chrono::high_resolution_clock::now();

			std::size_t num_min4 = generate_minimizers_17w(seq, 15, minimizers4.data());

			auto end = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration<double>(end - start).count();

			Assert::AreEqual(num_min1, num_min4);

			{
				std::stringstream message;
				message << "Vector version: " << duration << "s" << std::endl;
				Logger::WriteMessage(message.str().c_str());
			}

			//std::vector<uint32_t> minimizers2;
			//minimizers2.resize(4 * seq.size());

			//start = std::chrono::high_resolution_clock::now();

			//std::size_t num_min2 = generate_minimizers_q(seq, 15, 17, minimizers2.data());

			//end = std::chrono::high_resolution_clock::now();

			//duration = std::chrono::duration<double>(end - start).count();

			//Assert::AreEqual(num_min1, num_min2);

			//{
			//	std::stringstream message;
			//	message << "Deque version: " << duration << "s" << std::endl;
			//	Logger::WriteMessage(message.str().c_str());
			//}

			//for (std::size_t i = 0; i < num_min1; i++)
			//{
			//	Assert::AreEqual(minimizers2[i], minimizers4[i]);
			//}
		}
	};
}