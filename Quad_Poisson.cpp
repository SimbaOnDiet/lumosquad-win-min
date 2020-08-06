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
#include "QUAD_POISSON.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

QUAD_POISSON::QUAD_POISSON(int xRes, int yRes, int iterations) :
	_root(new CELL(1.0f, 1.0f, 0.0f, 0.0f))
{
	_root->refine();

	// figure out the max depth needed
	float xMax = log((float)xRes) / log(2.0f);
	float yMax = log((float)yRes) / log(2.0f);

	float max = (xMax > yMax) ? xMax : yMax;
	if (max - floor(max) > 1e-7)
		max = max + 1;
	max = floor(max);

	_maxRes = pow(2.0f, (float)max);
	_maxDepth = max;

	// create the blue noise
	_noiseFunc = new BLUE_NOISE(5.0f / (float)_maxRes);
	_noise = new bool[_maxRes * _maxRes];
	_noiseFunc->complete();
	_noiseFunc->maximize();
	_noiseFunc->writeToBool(_noise, _maxRes);

	_solver = new CG_SOLVER(_maxDepth, iterations);
}

QUAD_POISSON::~QUAD_POISSON()
{
	deleteGhosts();
	delete _root;
	delete _solver;
	delete _noiseFunc;
	delete[] _noise;
}


//////////////////////////////////////////////////////////////////////
// subdivide quadtree to max level for (xPos, yPos)
// returns a pointer to the cell that was created
//////////////////////////////////////////////////////////////////////
CELL* QUAD_POISSON::insert(float xPos, float yPos)
{
	int currentDepth = 0;
	CELL* currentCell = _root;
	bool existed = true;

	while (currentDepth < _maxDepth) {
		// find quadrant of current point
		float diff[2];
		diff[0] = xPos - currentCell->center[0];
		diff[1] = yPos - currentCell->center[1];
		int quadrant = 1;          //RA Note: I have to go back to look over it later...... It confused me. //RA Note: I have got it. But this style really confused me.
		if (diff[0] > 0.0f) {
			if (diff[1] < 0.0f)
				quadrant = 2;
		}
		else if (diff[1] < 0.0f)
			quadrant = 3;
		else
			quadrant = 0;

		// check if it exists
		if (currentCell->children[quadrant] == NULL) {
			existed = false;
			currentCell->refine();
		}

		// recurse to next level
		currentCell = currentCell->children[quadrant];

		// increment depth
		currentDepth++;
	}
	// if we had to subdivide to get the cell, add them to the list
	if (!existed)
		for (int i = 0; i < 4; i++)
		{
			_smallestLeaves.push_back(currentCell->parent->children[i]);
			setNoise(currentCell->parent->children[i]);
		}

	///////////////////////////////////////////////////////////////////
	// force orthogonal neighbors to be same depth
	// I have commented the first block, the rest follow the same flow
	///////////////////////////////////////////////////////////////////

	// see if neighbor exists
	CELL* north = currentCell->northNeighbor();
	if (north && north->depth != _maxDepth) {
		// while the neighbor needs to be refined
		while (north->depth != _maxDepth) {

			// refine it 
			north->refine();

			// set to the newly refined neighbor
			north = currentCell->northNeighbor();    //RA Note: A reliable line?
		}
		// add newly created nodes to the list
		for (int i = 0; i < 4; i++)
		{
			_smallestLeaves.push_back(north->parent->children[i]);
			setNoise(north->parent->children[i]);
		}
	}
	CELL* south = currentCell->southNeighbor();
	if (south && south->depth != _maxDepth) {
		while (south->depth != _maxDepth) {
			south->refine();
			south = currentCell->southNeighbor();
		}
		for (int i = 0; i < 4; i++)
		{
			_smallestLeaves.push_back(south->parent->children[i]);
			setNoise(south->parent->children[i]);
		}
	}
	CELL* west = currentCell->westNeighbor();
	if (west && west->depth != _maxDepth) {
		while (west->depth != _maxDepth) {
			west->refine();
			west = currentCell->westNeighbor();
		}
		for (int i = 0; i < 4; i++)
		{
			_smallestLeaves.push_back(west->parent->children[i]);
			setNoise(west->parent->children[i]);
		}
	}
	CELL* east = currentCell->eastNeighbor();
	if (east && east->depth != _maxDepth) {
		while (east->depth != _maxDepth) {
			east->refine();
			east = currentCell->eastNeighbor();
		}
		for (int i = 0; i < 4; i++)
		{
			_smallestLeaves.push_back(east->parent->children[i]);
			setNoise(east->parent->children[i]);
		}
	}

	///////////////////////////////////////////////////////////////////
	// force diagonal neighbors to be same depth
	// The same flow follows as above, except that it makes sure that
	// the 'north' and 'south' neighbors already exist
	///////////////////////////////////////////////////////////////////

	if (north) {
		CELL* northwest = north->westNeighbor();
		if (northwest && northwest->depth != _maxDepth) {
			while (northwest->depth != _maxDepth) {
				northwest->refine();
				northwest = northwest->children[2];       //RA Note: This seems a little more reasonable.
			}
			for (int i = 0; i < 4; i++)
			{
				_smallestLeaves.push_back(northwest->parent->children[i]);
				setNoise(northwest->parent->children[i]);
			}
		}
		CELL* northeast = north->eastNeighbor();
		if (northeast && northeast->depth != _maxDepth) {
			while (northeast->depth != _maxDepth) {
				northeast->refine();
				northeast = northeast->children[3];
			}
			for (int i = 0; i < 4; i++)
			{
				_smallestLeaves.push_back(northeast->parent->children[i]);
				setNoise(northeast->parent->children[i]);
			}
		}
	}
	if (south) {      //RA Note: Same routine.
		CELL* southwest = south->westNeighbor();
		if (southwest && southwest->depth != _maxDepth) {
			while (southwest->depth != _maxDepth) {
				southwest->refine();
				southwest = southwest->children[1];
			}
			for (int i = 0; i < 4; i++)
			{
				_smallestLeaves.push_back(southwest->parent->children[i]);
				setNoise(southwest->parent->children[i]);
			}
		}
		CELL* southeast = south->eastNeighbor();
		if (southeast && southeast->depth != _maxDepth) {
			while (southeast->depth != _maxDepth) {
				southeast->refine();
				southeast = southeast->children[0];
			}
			for (int i = 0; i < 4; i++)
			{
				_smallestLeaves.push_back(southeast->parent->children[i]);
				setNoise(southeast->parent->children[i]);
			}
		}
	}

	return currentCell;
}


//////////////////////////////////////////////////////////////////////
// insert all leaves into a list
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::getAllLeaves(list<CELL*>& leaves, CELL* currentCell)
{
	// if we're at the root
	if (currentCell == NULL)
	{
		getAllLeaves(leaves, _root);
		return;
	}

	// if we're at a leaf, add it to the list
	if (currentCell->children[0] == NULL)
	{
		leaves.push_back(currentCell);
		return;
	}

	// if children exist, call recursively
	for (int x = 0; x < 4; x++)
		getAllLeaves(leaves, currentCell->children[x]);
}

//////////////////////////////////////////////////////////////////////
// insert all leaves not on the boundary into a list
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::getEmptyLeaves(list<CELL*>& leaves, CELL* currentCell)
{
	// if we're at the root
	if (currentCell == NULL) {
		getEmptyLeaves(leaves, _root);
		return;
	}

	// if we're at a leaf, check if it's a boundary and then
	// add it to the list
	if (currentCell->children[0] == NULL) {
		if (!(currentCell->boundary))
			leaves.push_back(currentCell);
		return;
	}

	// if children exist, call recursively
	for (int x = 0; x < 4; x++)
		getEmptyLeaves(leaves, currentCell->children[x]);
}

//////////////////////////////////////////////////////////////////////
// balance the current tree
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::balance()
{
	// collect all the leaf nodes
	list<CELL*> leaves;
	getAllLeaves(leaves);

	// while the list is not empty
	list<CELL*>::iterator cellIterator = leaves.begin();
	for (cellIterator = leaves.begin(); cellIterator != leaves.end(); cellIterator++) {
		CELL* currentCell = *cellIterator;

		// if a north neighbor exists
		CELL* north = currentCell->northNeighbor();
		if (north != NULL)
			// while the neighbor is not balanced
			while (north->depth < currentCell->depth - 1) {
				// refine it
				north->refine();

				// add the newly refined nodes to the list of
				// those to be checked
				for (int x = 0; x < 4; x++)
					leaves.push_back(north->children[x]);

				// set the cell to the newly created one
				north = currentCell->northNeighbor();
			}

		// the rest of the blocks flow the same as above
		CELL* south = currentCell->southNeighbor();
		if (south != NULL)
			while (south->depth < currentCell->depth - 1) {
				south->refine();
				for (int x = 0; x < 4; x++)
					leaves.push_back(south->children[x]);
				south = currentCell->southNeighbor();
			}

		CELL* west = currentCell->westNeighbor();
		if (west != NULL)
			while (west->depth < currentCell->depth - 1) {
				west->refine();
				for (int x = 0; x < 4; x++)
					leaves.push_back(west->children[x]);
				west = currentCell->westNeighbor();
			}

		CELL* east = currentCell->eastNeighbor();
		if (east != NULL)
			while (east->depth < currentCell->depth - 1) {
				east->refine();
				for (int x = 0; x < 4; x++)
					leaves.push_back(east->children[x]);
				east = currentCell->eastNeighbor();
			}
	}
}

//////////////////////////////////////////////////////////////////////
// build the neighbor lists of the current quadtree
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::buildNeighbors()
{
	balance();

	// collect all the leaf nodes
	list<CELL*> leaves;
	getAllLeaves(leaves);

	list<CELL*>::iterator cellIterator = leaves.begin();
	for (cellIterator = leaves.begin(); cellIterator != leaves.end(); cellIterator++)
	{
		CELL* currentCell = *cellIterator;

		// build north neighbors
		CELL* north = currentCell->northNeighbor();
		if (north != NULL) {
			if (north->children[0] == NULL) {
				currentCell->neighbors[0] = north;
				currentCell->neighbors[1] = NULL;
			}
			else {
				currentCell->neighbors[0] = north->children[3];
				currentCell->neighbors[1] = north->children[2];
			}
		}
		// else build a ghost cell
		else
			currentCell->neighbors[0] = new CELL(currentCell->depth);

		// build east neighbors
		CELL* east = currentCell->eastNeighbor();
		if (east != NULL) {
			if (east->children[0] == NULL) {
				currentCell->neighbors[2] = east;
				currentCell->neighbors[3] = NULL;
			}
			else {
				currentCell->neighbors[2] = east->children[0];
				currentCell->neighbors[3] = east->children[3];
			}
		}
		// else build a ghost cell
		else
			currentCell->neighbors[2] = new CELL(currentCell->depth);

		// build south neighbors
		CELL* south = currentCell->southNeighbor();
		if (south != NULL) {
			if (south->children[0] == NULL) {
				currentCell->neighbors[4] = south;
				currentCell->neighbors[5] = NULL;
			}
			else {
				currentCell->neighbors[4] = south->children[1];
				currentCell->neighbors[5] = south->children[0];
			}
		}
		// else build a ghost cell
		else
			currentCell->neighbors[4] = new CELL(currentCell->depth);

		// build west neighbors
		CELL* west = currentCell->westNeighbor();
		if (west != NULL) {
			if (west->children[0] == NULL) {
				currentCell->neighbors[6] = west;
				currentCell->neighbors[7] = NULL;
			}
			else {
				currentCell->neighbors[6] = west->children[2];
				currentCell->neighbors[7] = west->children[1];
			}
		}
		// else build a ghost cell
		else
			currentCell->neighbors[6] = new CELL(currentCell->depth);
	}
}

//////////////////////////////////////////////////////////////////////
// delete ghost cells
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::deleteGhosts(CELL* currentCell)
{
	// if at the base, call on the root
	if (currentCell == NULL)
	{
		deleteGhosts(_root);
		return;
	}

	// if there are children, delete those too
	if (currentCell->children[0]) {
		// call recursively
		for (int x = 0; x < 4; x++)
			deleteGhosts(currentCell->children[x]);
		return;
	}

	// check the neighbors for stuff to delete
	#pragma omp parallel for
	for (int x = 0; x < 8; x++) {
		// if the neighbor exists
		if (currentCell->neighbors[x])
			// and if it is a ghost cell, delete it
			if (currentCell->neighbors[x]->parent == NULL)
				delete currentCell->neighbors[x];
	}
}

//////////////////////////////////////////////////////////////////////
// solve the Poisson problem
//////////////////////////////////////////////////////////////////////
int QUAD_POISSON::solve() {
	// maintain the quadtree
	balance();
	buildNeighbors();

	// retrieve leaves at the lowest level
	_emptyLeaves.clear();
	getEmptyLeaves(_emptyLeaves);

	static bool firstSolve = true;
	static int iterations = _solver->iterations();

	// do a full precision solve the first time
	if (firstSolve)   //RA Note: Check it up later.
	{
		iterations = _solver->iterations();
		_solver->iterations() = 10000;
		firstSolve = false;
	}
	else
		_solver->iterations() = iterations;

	// return the number of iterations
	return _solver->solve(_emptyLeaves);
};


//////////////////////////////////////////////////////////////////////
// get the leafnode that corresponds to the coordinate
//////////////////////////////////////////////////////////////////////
CELL* QUAD_POISSON::getLeaf(float xPos, float yPos)
{
	CELL* currentCell = _root;

	while (currentCell->children[0] != NULL)
	{
		// find quadrant of current point
		float diff[2];
		diff[0] = xPos - currentCell->center[0];
		diff[1] = yPos - currentCell->center[1];
		int quadrant = 1;
		if (diff[0] > 0.0f)
		{
			if (diff[1] < 0.0f)
				quadrant = 2;
		}
		else if (diff[1] < 0.0f)
			quadrant = 3;
		else
			quadrant = 0;

		// check if it exists
		if (currentCell->children[quadrant] != NULL)
			currentCell = currentCell->children[quadrant];
	}
	return currentCell;
}

//////////////////////////////////////////////////////////////////////
// check if a cell hits a noise node
//////////////////////////////////////////////////////////////////////
void QUAD_POISSON::setNoise(CELL* cell)
{
	if (!(cell->state == EMPTY))
		return;

	int x = cell->center[0] * _maxRes;
	int y = cell->center[1] * _maxRes;

	if (_noise[x + y * _maxRes])          //RA Note: Check it up later.
	{
		cell->boundary = true;
		cell->state = ATTRACTOR;
		cell->potential = 0.5f;
		cell->candidate = true;
	}
}
