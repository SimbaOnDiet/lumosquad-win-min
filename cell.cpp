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
#include "cell.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// normal cell constructor
CELL::CELL(float north, float east, float south, float west, CELL* parent, int depth) :
	parent(parent), depth(depth), index(-1), candidate(false),
	boundary(false), potential(0.0f), state(EMPTY)
{
	#pragma omp parallel for
	for (int x = 0; x < 4; x++)
		children[x] = NULL;

	#pragma omp parallel for
	for (int x = 0; x < 8; x++)
		neighbors[x] = NULL;

	bounds[0] = north; bounds[1] = east; bounds[2] = south; bounds[3] = west;

	center[0] = (bounds[1] + bounds[3]) * 0.5f;
	center[1] = (bounds[0] + bounds[2]) * 0.5f;
}

// ghost cell constructor
CELL::CELL(int depth) : parent(NULL), depth(depth), index(-1), candidate(false),
boundary(true), potential(0.0f), state(EMPTY)
{
	#pragma omp parallel for
	for (int x = 0; x < 4; x++)
		children[x] = NULL;

	#pragma omp parallel for
	for (int x = 0; x < 8; x++)
		neighbors[x] = NULL;

	bounds[0] = 0.0f; bounds[1] = 0.0f; bounds[2] = 0.0f; bounds[3] = 0.0f;
	center[0] = 0.0f; center[1] = 0.0f;
}

CELL::~CELL() {
	int x;

	for (x = 0; x < 4; x++)
		if (children[x] != NULL)
		{
			delete children[x];
			children[x] = NULL;
		}
}

//////////////////////////////////////////////////////////////////////
// refine current cell
//////////////////////////////////////////////////////////////////////
void CELL::refine() {
	if (children[0] != NULL) return;
	float center[] = { (bounds[0] + bounds[2]) * 0.5f, (bounds[1] + bounds[3]) * 0.5f };

	children[0] = new CELL(bounds[0], center[1], center[0], bounds[3], this, depth + 1);
	children[1] = new CELL(bounds[0], bounds[1], center[0], center[1], this, depth + 1);
	children[2] = new CELL(center[0], bounds[1], bounds[2], center[1], this, depth + 1);
	children[3] = new CELL(center[0], center[1], bounds[2], bounds[3], this, depth + 1);

	children[0]->potential = potential;
	children[1]->potential = potential;
	children[2]->potential = potential;
	children[3]->potential = potential;
}

//////////////////////////////////////////////////////////////////////
// return north neighbor to current cell
//////////////////////////////////////////////////////////////////////
CELL* CELL::northNeighbor()
{
	// if it is the root
	if (this->parent == NULL) return NULL;

	// if it is the southern child of the parent
	if (parent->children[3] == this) return parent->children[0];
	if (parent->children[2] == this) return parent->children[1];

	// else look up higher
	CELL* mu = parent->northNeighbor();

	// if there are no more children to look at,
	// this is the answer
	if (mu == NULL || mu->children[0] == NULL) return mu;
	// if it is the NW child of the parent
	else if (parent->children[0] == this) return mu->children[3];
	// if it is the NE child of the parent
	else return mu->children[2];
}

//////////////////////////////////////////////////////////////////////
// return north neighbor to current cell
//////////////////////////////////////////////////////////////////////
CELL* CELL::southNeighbor()
{
	// if it is the root
	if (this->parent == NULL) return NULL;

	// if it is the northern child of the parent
	if (parent->children[0] == this) return parent->children[3];
	if (parent->children[1] == this) return parent->children[2];

	// else look up higher
	CELL* mu = parent->southNeighbor();

	// if there are no more children to look at,
	// this is the answer
	if (mu == NULL || mu->children[0] == NULL) return mu;
	// if it is the SW child of the parent
	else if (parent->children[3] == this) return mu->children[0];
	// if it is the SE child of the parent
	else return mu->children[1];
}

//////////////////////////////////////////////////////////////////////
// return north neighbor to current cell
//////////////////////////////////////////////////////////////////////
CELL* CELL::westNeighbor()
{
	// if it is the root
	if (this->parent == NULL) return NULL;

	// if it is the eastern child of the parent
	if (parent->children[1] == this) return parent->children[0];
	if (parent->children[2] == this) return parent->children[3];

	// else look up higher
	CELL* mu = parent->westNeighbor();

	// if there are no more children to look at,
	// this is the answer
	if (mu == NULL || mu->children[0] == NULL) return mu;
	// if it is the NW child of the parent
	else if (parent->children[0] == this) return mu->children[1];
	// if it is the SW child of the parent
	else return mu->children[2];
}

//////////////////////////////////////////////////////////////////////
// return north neighbor to current cell
//////////////////////////////////////////////////////////////////////
CELL* CELL::eastNeighbor()
{
	// if it is the root
	if (this->parent == NULL) return NULL;

	// if it is the western child of the parent
	if (parent->children[0] == this) return parent->children[1];
	if (parent->children[3] == this) return parent->children[2];

	// else look up higher
	CELL* mu = parent->eastNeighbor();

	// if there are no more children to look at,
	// this is the answer
	if (mu == NULL || mu->children[0] == NULL) return mu;
	// if it is the NE child of the parent
	else if (parent->children[1] == this) return mu->children[0];
	// if it is the SE child of the parent
	else return mu->children[3];
}

