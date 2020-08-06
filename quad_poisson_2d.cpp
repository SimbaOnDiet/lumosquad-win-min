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
#include "FFT.h"
#include "APSF.h"
#include "QUAD_poisson_2D.h"

////////////////////////////////////////////////////////////////////////////
// render the glow
////////////////////////////////////////////////////////////////////////////
void QUAD_DBM_2D::renderglow(void)
{
	const int scale = 1;
	int w = _xRes * scale;
	int h = _yRes * scale;

	// draw the DAG
	float*& source = renderOffscreen(scale);

	// if there is no input dimensions specified, else there were input
	// image dimensions, so crop it
	
	// copy out the cropped version
	int wCropped = w * scale;
	int hCropped = h * scale;
	float* cropped = new float[wCropped * hCropped];
	cout << endl << " Generating rendered image width: " << wCropped << " height: " << hCropped << endl;

	#pragma omp parallel for
	for (int y = 0; y < hCropped; y++)
		for (int x = 0; x < wCropped; x++)
		{
			int uncroppedIndex = x + y * w;
			int croppedIndex = x + y * wCropped;
			cropped[croppedIndex] = source[uncroppedIndex];
		}

	// create the filter
	APSF apsf(w);
	apsf.generateKernelFast();

	// convolve with FFT
	bool success = FFT::convolve(cropped, apsf.kernel(), wCropped, hCropped, apsf.res(), apsf.res());

	
	if (success) {
		Alltofile(cropped,w);
	}
	else
		cout << " Final image generation failed." << endl;
	


	delete[] cropped;
}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

QUAD_DBM_2D::QUAD_DBM_2D(int xRes, int yRes, int iterations) :
	_xRes(xRes),
	_yRes(yRes),
	_bottomHit(0),
	_iterations(iterations),
	_quadPoisson(NULL),
	_skips(10),
	_twister(123456)
{
	allocate();
	_dag = new DAG(_xRes, _yRes);

	// calculate dimensions
	_dx = 1.0f / (float)_xRes;
	_dy = 1.0f / (float)_yRes;
	if (_dx < _dy) _dy = _dx;
	else _dx = _dy;

	_maxRes = _xRes * _yRes;
}

QUAD_DBM_2D::~QUAD_DBM_2D()
{
	deallocate();
}

void QUAD_DBM_2D::allocate()
{
	_quadPoisson = new QUAD_POISSON(_xRes, _yRes, _iterations);
	_xRes = _yRes = _quadPoisson->maxRes();
}

void QUAD_DBM_2D::deallocate()
{
	if (_quadPoisson) delete _quadPoisson;
}

//////////////////////////////////////////////////////////////////////
// check neighbors for any candidate nodes
//////////////////////////////////////////////////////////////////////
void QUAD_DBM_2D::checkForCandidates(CELL* cell)
{
	int maxDepth = _quadPoisson->maxDepth();

	CELL* north = cell->northNeighbor();
	if (north) {
		if (north->depth == maxDepth) {
			if (!north->candidate) {
				_candidates.push_back(north);
				north->candidate = true;
			}

			CELL* northeast = north->eastNeighbor();
			if (northeast && !northeast->candidate) {
				_candidates.push_back(northeast);
				northeast->candidate = true;
			}
			CELL* northwest = north->westNeighbor();
			if (northwest && !northwest->candidate) {
				_candidates.push_back(northwest);
				northwest->candidate = true;
			}
		}
	}

	CELL* east = cell->eastNeighbor();
	if (east && !east->candidate) {
		_candidates.push_back(east);
		east->candidate = true;
	}

	CELL* south = cell->southNeighbor();
	if (south) {
		if (!south->candidate) {
			_candidates.push_back(south);
			south->candidate = true;
		}

		CELL* southeast = south->eastNeighbor();
		if (southeast && !southeast->candidate) {
			_candidates.push_back(southeast);
			southeast->candidate = true;
		}

		CELL* southwest = south->westNeighbor();
		if (southwest && !southwest->candidate) {
			_candidates.push_back(southwest);
			southwest->candidate = true;
		}
	}

	CELL* west = cell->westNeighbor();
	if (west && !west->candidate) {
		_candidates.push_back(west);
		west->candidate = true;
	}
}

//////////////////////////////////////////////////////////////////////
// add particle to the aggregate
//////////////////////////////////////////////////////////////////////
bool QUAD_DBM_2D::addParticle()
{
	static float invSqrtTwo = 1.0f / sqrt(2.0f);
	static int totalParticles = 0;
	static int skipSolve = 0;

	// compute the potential
	int iterations = 0;
	if (!skipSolve)
		iterations = _quadPoisson->solve();
	skipSolve++;
	if (skipSolve == _skips) skipSolve = 0;

	// construct probability distribution
	vector<float> probabilities;
	float totalPotential = 0.0f;
	for (int x = 0; x < _candidates.size(); x++)
	{
		if (_candidates[x]->candidate)
		{
			probabilities.push_back(_candidates[x]->potential);
			totalPotential += _candidates[x]->potential;
		}
		else
			probabilities.push_back(0.0f);
	}

	// get all the candidates
	// if none are left, stop
	if (_candidates.size() == 0) {
		return false;
	}

	// if there is not enough potential, go Brownian
	int toAddIndex = 0;
	if (totalPotential < 1e-8)
		toAddIndex = _candidates.size() * _twister.getDoubleLR();
	// else follow DBM algorithm
	else
	{
		// add a neighbor
		float random = _twister.getDoubleLR();
		float invTotalPotential = 1.0f / totalPotential;
		float potentialSeen = probabilities[0] * invTotalPotential;
		while ((potentialSeen < random) && (toAddIndex < _candidates.size()))
		{
			toAddIndex++;
			potentialSeen += probabilities[toAddIndex] * invTotalPotential;
		}
	}
	_candidates[toAddIndex]->boundary = true;
	_candidates[toAddIndex]->potential = 0.0f;
	_candidates[toAddIndex]->state = NEGATIVE;

	CELL* neighbor = NULL;
	CELL* added = _candidates[toAddIndex];
	CELL* north = added->northNeighbor();
	if (north)
	{
		if (north->state == NEGATIVE)
			neighbor = north;
		CELL* northeast = north->eastNeighbor();
		if (northeast && northeast->state == NEGATIVE)
			neighbor = northeast;
		CELL* northwest = north->westNeighbor();
		if (northwest && northwest->state == NEGATIVE)
			neighbor = northwest;
	}
	CELL* east = added->eastNeighbor();
	if (east && east->state == NEGATIVE)
		neighbor = east;

	CELL* south = added->southNeighbor();
	if (south)
	{
		if (south->state == NEGATIVE)
			neighbor = south;
		CELL* southeast = south->eastNeighbor();
		if (southeast && southeast->state == NEGATIVE)
			neighbor = southeast;
		CELL* southwest = south->westNeighbor();
		if (southwest && southwest->state == NEGATIVE)
			neighbor = southwest;
	}
	CELL* west = added->westNeighbor();
	if (west && west->state == NEGATIVE)
		neighbor = west;

	// insert it as a node for bookkeeping
	_quadPoisson->insert(added->center[0], added->center[1]);
	checkForCandidates(added);

	// insert into the DAG
	int newIndex = (int)(added->center[0] * _xRes) +
		(int)(added->center[1] * _yRes) * _xRes;
	int neighborIndex = (int)(neighbor->center[0] * _xRes) +
		(int)(neighbor->center[1] * _yRes) * _xRes;
	_dag->addSegment(newIndex, neighborIndex);

	totalParticles++;
	if (!(totalParticles % 200))
		cout << " " << totalParticles;

	//_dag->drawNode();

	//Little change here
	bool stillgoing = true;
	stillgoing=hitGround(added);

	return !stillgoing;
}

//////////////////////////////////////////////////////////////////////
// hit ground yet?
//////////////////////////////////////////////////////////////////////
bool QUAD_DBM_2D::hitGround(CELL* cell)
{
	if (_bottomHit)
		return true;

	if (!cell)
		return false;

	bool hit = false;
	if (cell->northNeighbor())
	{
		CELL* north = cell->northNeighbor();
		if (north->state == POSITIVE)
			hit = true;
		if (north->eastNeighbor()->state == POSITIVE)
			hit = true;
		if (north->westNeighbor()->state == POSITIVE)
			hit = true;
	}
	if (cell->eastNeighbor())
		if (cell->eastNeighbor()->state == POSITIVE)
			hit = true;
	if (cell->southNeighbor())
	{
		CELL* south = cell->southNeighbor();
		if (south->state == POSITIVE)
			hit = true;
		if (south->eastNeighbor()->state == POSITIVE)
			hit = true;
		if (south->westNeighbor()->state == POSITIVE)
			hit = true;
	}
	if (cell->westNeighbor())
		if (cell->westNeighbor()->state == POSITIVE)
			hit = true;

	if (hit)
	{
		cout << endl<<"Hit Ground!"<<endl;
		_bottomHit = (int)(cell->center[0] * _xRes) +
			(int)(cell->center[1] * _yRes) * _xRes;
		//_dag = new DAG(_xRes, _yRes);
		_dag->buildLeader(_bottomHit);
		renderglow();
		return true;
	}
	return false;
}

void QUAD_DBM_2D::Alltofile(float* points_, int Res) {
	string filename = "rendered";
	ofstream outfile(filename);

	for (int i = 0; i < Res; i++) {
		for (int j = 0; j < Res; j++) {
			outfile << points_[i*Res + j] << " ";
		}
		outfile << endl;
	}
	outfile.close();


	delete[] points_;

	cout << "Completed! You can close this program now." << endl;
}

////////////////////////////////////////////////////////////////////
// read in attractors from an image
////////////////////////////////////////////////////////////////////
bool QUAD_DBM_2D::readImage(unsigned char* initial,
	unsigned char* attractors,
	unsigned char* repulsors,
	unsigned char* terminators,
	int xRes, int yRes)
{


	bool initialFound = false;
	bool terminateFound = false;

	int index = 0;
	for (int y = 0; y < yRes; y++)
		for (int x = 0; x < xRes; x++, index++)
		{
			// insert initial condition
			if (initial[index])
			{
				// insert something
				CELL* negative = _quadPoisson->insert(x, y);
				negative->boundary = true;
				negative->potential = 0.0f;
				negative->state = NEGATIVE;
				negative->candidate = true;

				checkForCandidates(negative);

				initialFound = true;
			}

			// insert attractors
			if (attractors[index])
			{
				// insert something
				CELL* positive = _quadPoisson->insert(x, y);
				positive->boundary = true;
				positive->potential = 1.0f;
				positive->state = ATTRACTOR;
				positive->candidate = true;
			}

			// insert repulsors
			if (repulsors[index])
			{
				// only insert the repulsor if it is the edge of a repulsor
				bool edge = false;

				if (x != 0)
				{
					if (!repulsors[index - 1]) edge = true;
					if (y != 0 && !repulsors[index - xRes - 1]) edge = true;
					if (y != yRes - 1 && !repulsors[index + xRes - 1]) edge = true;
				}
				if (x != _xRes - 1)
				{
					if (!repulsors[index + 1]) edge = true;
					if (y != 0 && !repulsors[index - xRes + 1]) edge = true;
					if (y != yRes - 1 && !repulsors[index + xRes + 1]) edge = true;
				}
				if (y != 0 && !repulsors[index - xRes]) edge = true;
				if (y != yRes - 1 && !repulsors[index + xRes]) edge = true;

				if (edge)
				{
					// insert something
					CELL* negative = _quadPoisson->insert(x, y);
					negative->boundary = true;
					negative->potential = 0.0f;
					negative->state = REPULSOR;
					negative->candidate = true;
				}
			}

			// insert terminators
			if (terminators[index])
			{
				// insert something
				CELL* positive = _quadPoisson->insert(x, y);
				positive->boundary = true;
				positive->potential = 1.0f;
				positive->state = POSITIVE;
				positive->candidate = true;

				terminateFound = true;
			}
		}

	if (!initialFound) {
		cout << " The lightning does not start anywhere! " << endl;
		return false;
	}
	if (!terminateFound) {
		cout << " The lightning does not end anywhere! " << endl;
		return false;
	}

	return true;
}
