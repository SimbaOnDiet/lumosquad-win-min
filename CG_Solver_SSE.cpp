//////////////////////////////////////////////////////////////////////////////////
//
// THOR - A Phyciscal Based Lightning Generator
//
// A Program for the (Undergraduate Course) Electromagnetism Presentation
//
// Modified from LumosQuad by 
// Liang Yien, 
// Student,
// University of Science and Technology of China
// 
///////////////////////////////////////////////////////////////////////////////////
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  I makes no representations about the suitability of this software for
//  any purpose. It is provided "as is" without express or implied warranty.
//  
//  Permission to use, copy, modify and distribute this software and its
//  documentation for educational, research and non-profit purposes, without
//  fee, and without a written agreement is hereby granted, provided that the
//  above copyright notice and the following three paragraphs appear in all
//  copies.
//
//  I SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//   THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
//  PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN  "AS IS" BASIS, AND 
//  I HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS,
//  OR MODIFICATIONS.
//
///////////////////////////////////////////////////////////////////////////////////
//
//  This program is largely based on LumosQuad, which has the following copyright 
//  notice and three other free software copyright notice which LumosQuad uses.
//  However, for the reason that I would only like to keep the basic function of
//  generating lightning and make it easy to make others understand the source of 
//  this program, some of them are not included in this program. But I will just
//  keep their copyright notice anyway, just to thanks for their contribution to
//  the original LumosQuad.
//
//  I express my sincere  gratitude to Theodore Kim and Ming C. Lin, the author 
//  of LumosQuad,a wonderful software, and Fast Animation of Lightning Using an 
//  Adaptive Mesh, an inspiring paper, without which I could not have complete 
//  this program.
//  
//  This program also use Mersenne Twister Algorithm to generate random numbers,
//  which is attributed to Makoto Matsumoto and Takuji Nishimura. Its copyright 
//  notice is shown below.
//  
//  
///////////////////////////////////////////////////////////////////////////////////
//
// LumosQuad - A Lightning Generator
// Copyright 2007
// The University of North Carolina at Chapel Hill
//
//  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY WARRANTIES,
//  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
//  FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN
//  "AS IS" BASIS, AND THE UNIVERSITY OF NORTH CAROLINA HAS NO OBLIGATION TO
//  PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
//
//  Please send questions and comments about LumosQuad to kim@cs.unc.edu.
// 
///////////////////////////////////////////////////////////////////////////////////
//  OpenEXR Copyright Notice:
//
//  Copyright (c) 2002, Industrial Light & Magic, a division of Lucas
//  Digital Ltd. LLC
// 
//  All rights reserved.
// 
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//  *       Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  *       Redistributions in binary form must reproduce the above
//  copyright notice, this list of conditions and the following disclaimer
//  in the documentation and/or other materials provided with the
//  distribution.
//  *       Neither the name of Industrial Light & Magic nor the names of
//  its contributors may be used to endorse or promote products derived
//  from this software without specific prior written permission. 
// 
////////////////////////////////////////////////////////////////////////////////////
//
// This program contains a very thin wrapper to Daniel Dunbar's blue noise generator.
//
// For the original, untainted code, see: 
//   http://www.cs.virginia.edu/~gfx/pubs/antimony/
//
///////////////////////////////////////////////////////////////////////////////////
// This program used Mersenne Twister Random Number Generator, witten by 
// Takuji Nishimura and Makoto Matsumoto, and modified to be a C++ class 
// by Daniel Dunbar. Its copyright notice is shown below:
//
//   A C-program for MT19937, with initialization improved 2002/1/26.
//   Coded by Takuji Nishimura and Makoto Matsumoto.
//   Modified to be a C++ class by Daniel Dunbar.
//
//   Before using, initialize the state by using init_genrand(seed)  
//   or init_by_array(init_key, key_length).
//
//  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
//  All rights reserved.                          
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     3. The names of its contributors may not be used to endorse or promote 
//        products derived from this software without specific prior written 
//        permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
//   Any feedback is very welcome.
//   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
//
//////////////////////////////////////////////////////////////////////////
//  This program also uses FFTW. Here`s its copyright notice:
//
//
// Copyright (c) 2003, 2007-14 Matteo Frigo
// Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// FFTW, online at http://www.fftw.org
//
////////////////////////////////////////////////////////////////////////////




#include "stdafx.h"
#include "CG_Solver_SSE.h"
#include <stdlib.h>
#include <cstdlib>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CG_SOLVER_SSE::CG_SOLVER_SSE(int maxDepth, int iterations, int digits) :
	CG_SOLVER(maxDepth, iterations, digits)
{
}

CG_SOLVER_SSE::~CG_SOLVER_SSE()
{
	if (_direction) _aligned_free(_direction);
	if (_potential) _aligned_free(_potential);
	if (_residual)  _aligned_free(_residual);
	if (_q)         _aligned_free(_q);

	_direction = NULL;
	_residual = NULL;
	_q = NULL;
	_potential = NULL;
}

//////////////////////////////////////////////////////////////////////
// reallocate the sse arrays
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::reallocate()
{
	if (_arraySize >= _listSize) return;
	_arraySize = _listSize * 2;

	if (_arraySize % 4)
		_arraySize += 4 - _arraySize % 4;

	if (_direction) _aligned_free(_direction);
	if (_potential) _aligned_free(_potential);
	if (_residual)  _aligned_free(_residual);
	if (_q)         _aligned_free(_q);

	_direction = (float*)_aligned_malloc(_arraySize * sizeof(float), 16);
	_potential = (float*)_aligned_malloc(_arraySize * sizeof(float), 16);
	_residual = (float*)_aligned_malloc(_arraySize * sizeof(float), 16);
	_q = (float*)_aligned_malloc(_arraySize * sizeof(float), 16);

	return;
}

//////////////////////////////////////////////////////////////////////
// solve the linear system
//////////////////////////////////////////////////////////////////////
int CG_SOLVER_SSE::solve(list<CELL*> cells)
{
	// counters
	int x, y, index;
	list<CELL*>::iterator cellIterator;

	// i = 0
	int i = 0;

	// precalculate stencils
	calcStencils(cells);

	// reallocate scratch arrays if necessary
	_listSize = cells.size();
	reallocate();
	wipeSSE(_potential);
	wipeSSE(_direction);
	wipeSSE(_residual);
	wipeSSE(_q);

	// compute a new lexicographical order
	cellIterator = cells.begin();

	#pragma omp parallel for
	for (x = 0; x < _listSize; x++, cellIterator++)
	{
		CELL* cell = *cellIterator;
		cell->index = x;
		_potential[x] = cell->potential;
	}

	// r = b - Ax
	calcResidual(cells);

	// d = r
	copySSE(_direction, _residual);

	// deltaNew = r^T r
	float deltaNew = dotSSE(_residual, _residual);

	// delta0 = deltaNew
	float delta0 = deltaNew;

	// While deltaNew > (eps^2) * delta0
	float eps = pow(10.0f, (float)-_digits);
	float maxR = 2.0f * eps;
	while ((i < _iterations) && (maxR > eps))
	{
		// q = Ad
		cellIterator = cells.begin();
		#pragma omp parallel for
		for (y = 0; y < _listSize; y++, cellIterator++)
		{
			CELL* currentCell = *cellIterator;
			CELL** neighbors = currentCell->neighbors;
			float* stencil = currentCell->stencil;

			float neighborSum = 0.0f;
			for (int x = 0; x < 8; x++)
			{
				if (neighbors[x])
					neighborSum += _direction[neighbors[x]->index] * stencil[x];
			}
			_q[y] = -neighborSum + _direction[y] * currentCell->stencil[8];
		}

		// alpha = deltaNew / (transpose(d) * q)
		float alpha = dotSSE(_q, _direction);
		if (fabs(alpha) > 0.0f)
			alpha = deltaNew / alpha;

		// x = x + alpha * d
		saxpySSE(alpha, _direction, _potential);

		// r = r - alpha * q
		saxpySSE(-alpha, _q, _residual);
		maxR = maxSSE(_residual);

		// deltaOld = deltaNew
		float deltaOld = deltaNew;

		// deltaNew = transpose(r) * r
		deltaNew = dotSSE(_residual, _residual);

		// beta = deltaNew / deltaOld
		float beta = deltaNew / deltaOld;

		// d = r + beta * d
		saypxSSE(beta, _residual, _direction);

		// i = i + 1
		i++;
	}

	// copy back into the tree
	cellIterator = cells.begin();
	#pragma omp parallel for
	for (x = 0; x < _listSize; x++, cellIterator++)
		(*cellIterator)->potential = _potential[x];

	return i;
}

//////////////////////////////////////////////////////////////////////
// dot product of two vectors
//////////////////////////////////////////////////////////////////////
float CG_SOLVER_SSE::dotSSE(float* x, float* y)
{
	__m128 sum = _mm_set_ps1(0.0f);
	__m128* xSSE = (__m128*)x;
	__m128* ySSE = (__m128*)y;
	__m128 temp;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		temp = _mm_mul_ps(*xSSE, *ySSE);
		sum = _mm_add_ps(sum, temp);
		xSSE++;
		ySSE++;
	}

	union u {
		__m128 m;
		float f[4];
	} extract;
	extract.m = sum;
	return extract.f[0] + extract.f[1] + extract.f[2] + extract.f[3];
}

//////////////////////////////////////////////////////////////////////
// scalar 'a' x + y
// Y = aX + Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::saxpySSE(float s, float* x, float* y)
{
	__m128* ySSE = (__m128*)y;
	__m128* xSSE = (__m128*)x;
	__m128 sSSE = _mm_set_ps1(s);
	__m128 temp;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		temp = _mm_mul_ps(*xSSE, sSSE);
		*ySSE = _mm_add_ps(*ySSE, temp);

		xSSE++;
		ySSE++;
	}

}

//////////////////////////////////////////////////////////////////////
// scalar 'a' y + x
// Y = aY + X
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::saypxSSE(float s, float* x, float* y)
{
	__m128* ySSE = (__m128*)y;
	__m128* xSSE = (__m128*)x;
	__m128 sSSE = _mm_set_ps1(s);
	__m128 temp;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		temp = _mm_mul_ps(*ySSE, sSSE);
		*ySSE = _mm_add_ps(*xSSE, temp);

		xSSE++;
		ySSE++;
	}
}

//////////////////////////////////////////////////////////////////////
// scalar 'a' y + x
// Y = aY + X
//////////////////////////////////////////////////////////////////////
float CG_SOLVER_SSE::maxSSE(float* x)
{
	__m128 maxFoundSSE = _mm_set_ps1(0.0f);
	__m128* xSSE = (__m128*)x;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		maxFoundSSE = _mm_max_ps(*xSSE, maxFoundSSE);
		xSSE++;
	}

	union u {
		__m128 m;
		float f[4];
	} extract;
	extract.m = maxFoundSSE;
	float maxFound = extract.f[0] > extract.f[1] ? extract.f[0] : extract.f[1];
	maxFound = maxFound > extract.f[2] ? maxFound : extract.f[2];
	return maxFound > extract.f[3] ? maxFound : extract.f[3];
}

//////////////////////////////////////////////////////////////////////
// SSE add
// Y = X + Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::addSSE(float* x, float* y)
{
	__m128* ySSE = (__m128*)y;
	__m128* xSSE = (__m128*)x;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		*ySSE = _mm_add_ps(*ySSE, *xSSE);
		xSSE++;
		ySSE++;
	}
}

//////////////////////////////////////////////////////////////////////
// SSE multiply
// Y = X * Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::multiplySSE(float* x, float* y)
{
	__m128* ySSE = (__m128*)y;
	__m128* xSSE = (__m128*)x;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		*ySSE = _mm_mul_ps(*ySSE, *xSSE);
		xSSE++;
		ySSE++;
	}
}

//////////////////////////////////////////////////////////////////////
// SSE multiply
// Z = X * Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::multiplySSE(float* x, float* y, float* z)
{
	__m128* zSSE = (__m128*)z;
	__m128* ySSE = (__m128*)y;
	__m128* xSSE = (__m128*)x;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		*zSSE = _mm_mul_ps(*ySSE, *xSSE);
		xSSE++;
		ySSE++;
		zSSE++;
	}
}

//////////////////////////////////////////////////////////////////////
// SSE multiply
// Z = W - X * Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::multiplySubtractSSE(float* w, float* x, float* y, float* z)
{
	__m128* zSSE = (__m128*)z;
	__m128* ySSE = (__m128*)y;
	__m128* xSSE = (__m128*)x;
	__m128* wSSE = (__m128*)w;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		*zSSE = _mm_mul_ps(*ySSE, *xSSE);
		*zSSE = _mm_sub_ps(*wSSE, *zSSE);

		xSSE++;
		ySSE++;
		zSSE++;
		wSSE++;
	}
}

//////////////////////////////////////////////////////////////////////
// SSE set
// X = val
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::setSSE(float* x, float val)
{
	__m128* xSSE = (__m128*)x;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		*xSSE = _mm_set_ps1(val);
		xSSE++;
	}
}

//////////////////////////////////////////////////////////////////////
// SSE set
// X = 0
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::wipeSSE(float* x)
{
	__m128* xSSE = (__m128*)x;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		*xSSE = _mm_setzero_ps();
		xSSE++;
	}
}

//////////////////////////////////////////////////////////////////////
// SSE set
// X = Y
//////////////////////////////////////////////////////////////////////
void CG_SOLVER_SSE::copySSE(float* x, float* y)
{
	__m128* ySSE = (__m128*)y;

	#pragma omp parallel for
	for (int index = 0; index < _arraySize / 4; index++)
	{
		_mm_store_ps(x, *ySSE);
		x += 4;
		ySSE++;
	}
}
