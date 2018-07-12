/*
 * Original work Copyright (c) 2016, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COMPARATORS_HEADER_H
#define COMPARATORS_HEADER_H


#include <algorithm>
/**
 * Comparator used to order a list of patches and distances (small to big)
 */
class ComparePairs {
	public:
		bool operator()(std::pair<float,unsigned>& a, std::pair<float,unsigned>& b)
		{
			return a.first < b.first;
		}
};

/**
 * Comparator used to order a list of patches and distances (big to small)
 */
class ComparePairsInv {
	public:
		bool operator()(std::pair<float,unsigned>& a, std::pair<float,unsigned>& b)
		{
			return a.first > b.first;
		}
};

#endif
