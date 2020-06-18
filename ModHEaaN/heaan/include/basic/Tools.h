/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_BASIC_TOOLS_H_
#define HEAAN_BASIC_TOOLS_H_

#include "Define.h"

static void __heaan_bit_reverse_array(vector<word64> &vec)
{
	size_t size = vec.size();
	for (size_t i = 1, j = 0; i < size; ++i)
	{
		size_t bit = size >> 1;
		for (; j >= bit; bit >>= 1)
		{
			j -= bit;
		}
		j += bit;
		if (i < j)
		{
			swap(vec[i], vec[j]);
		}
	}
}

static void __heaan_bit_reverse_array(Message& vec)
{
	size_t size = vec.size();
	for (size_t i = 1, j = 0; i < size; ++i)
	{
		size_t bit = size >> 1;
		for (; j >= bit; bit >>= 1)
		{
			j -= bit;
		}
		j += bit;
		if (i < j)
		{
			swap(vec[i], vec[j]);
		}
	}
}

static void __heaan_find_prime_factors(vector<word64> &s, word64 number)
{
	while (number % 2 == 0)
	{
		s.push_back(2);
		number /= 2;
	}
	for (word64 i = 3; i < sqrt(number); i++)
	{
		while (number % i == 0)
		{
			s.push_back(i);
			number /= i;
		}
	}
	if (number > 2)
	{
		s.push_back(number);
	}
}

#endif
