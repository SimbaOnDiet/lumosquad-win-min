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
#include "dag.h"
#include "param.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DAG::DAG(int xRes, int yRes) :
	_xRes(xRes),
	_yRes(yRes),
	_width(xRes),
	_height(yRes),
	_dx(1.0f / xRes),
	_dy(1.0f / yRes),
	_root(NULL),
	_totalNodes(0),
	_bottomHit(-1),
	points(NULL),
	_secondaryIntensity(0.3f),
	_leaderIntensity(0.75f)
{
	(_dx < _dy) ? _dy = _dx : _dy = _dx;
	_offscreenBuffer = NULL;
}

DAG::~DAG()
{
	if (_offscreenBuffer) delete[] _offscreenBuffer;

	deleteNode(_root);
}

//////////////////////////////////////////////////////////////////////
// delete the nodes
//////////////////////////////////////////////////////////////////////
void DAG::deleteNode(NODE* root)
{
	if (root == NULL) return;

	for (unsigned int x = 0; x < root->neighbors.size(); x++)
		deleteNode(root->neighbors[x]);

	delete root;
	root = NULL;
}

//////////////////////////////////////////////////////////////////////
// add line segment 'index' to segment list
//////////////////////////////////////////////////////////////////////
bool DAG::addSegment(int index, int neighbor)
{
	if (_root)
	{
		// find corresponding root in DAG
		map<int, NODE*>::iterator iter = _hash.find(neighbor);
		NODE* root = (*iter).second;

		if (root == NULL) return false;

		// add to DAG
		NODE* newNode = new NODE(index);
		newNode->parent = root;
		root->neighbors.push_back(newNode);
		drawNode(newNode);

		// add to hash table
		_hash.insert(map<int, NODE*>::value_type(index, newNode));
	}
	else
	{
		// make the root
		_root = new NODE(neighbor);
		_hash.insert(map<int, NODE*>::value_type(neighbor, _root));

		// then do the add
		NODE* newNode = new NODE(index);
		newNode->parent = _root;
		_root->neighbors.push_back(newNode);
		drawNode(newNode);

		// add to hash table
		_hash.insert(map<int, NODE*>::value_type(index, newNode));
	}

	_totalNodes++;
	/*
	if (_totalNodes == 1) {
		cout << "  [Now calls AddSergment]  " << endl;
	}
	if (!(_totalNodes % 200)) {
		cout << "(Sergent)" << _totalNodes << " ";
	}
	*/

	return true;
}

//////////////////////////////////////////////////////////////////////
// build the leader chain
//////////////////////////////////////////////////////////////////////
void DAG::buildLeader(int bottomHit)
{
	_bottomHit = bottomHit;

	// get pointer to bottommost node
	map<int, NODE*>::iterator iter = _hash.find(bottomHit);
	NODE* child = (*iter).second;

	// crawl up the tree
	if (child == NULL) cout << "WRONG!" << endl;

	while (child != NULL)
	{
		// tag segment
		child->leader = true;
		child->secondary = false;

		// look for side branches
		for (unsigned int x = 0; x < child->neighbors.size(); x++)
			if (!(child->neighbors[x]->leader))
				buildBranch(child->neighbors[x], 1);

		// advance child
		child = child->parent;
	}
	cout <<  "Root found!" << endl;
	buildIntensity(_root);
	int Res_ = xRes();   //BE CAREFUL HERE
	float* data= drawAll(_root, Res_, NULL);  //Res
	Alltofile(data,Res_);
}

//////////////////////////////////////////////////////////////////////
// build the side branch
//////////////////////////////////////////////////////////////////////
void DAG::buildBranch(NODE* node, int depth)
{
	node->depth = depth;
	node->leader = false;

	static int count = 0;
	count++;
	if (count == 1) cout << "Now building branches:" << endl;
	if (!(count % 200)) cout << count << "  ";

	// look for side branches
	for (unsigned int x = 0; x < node->neighbors.size(); x++)
		if (!(node->neighbors[x]->leader))
			buildBranch(node->neighbors[x], depth + 1);
}

//////////////////////////////////////////////////////////////////////
// DAG the DAG segments
//////////////////////////////////////////////////////////////////////
void DAG::drawNode(NODE* node)
{
	if (node == NULL) return;

	// DAG segments
	//int beginIndex = node->index;
	//int begin[2];
	//int end[2];
	//float dWidth = _dx;
	//float dHeight = _dy;

	//for (int x = 0; x < root->neighbors.size(); x++)
	//{
		// get end node
		//NODE* endNode = root->neighbors[x];
		//int endIndex = endNode->index;

		// DAG segments
		//begin[1] = beginIndex / _xRes;
		//begin[0] = (beginIndex - begin[1] * _xRes)/_xRes;
		//end[1] = endIndex / _xRes;
		//end[0] = endIndex - end[1] * _xRes;

		if (_bottomHit != -1)
			delete[] refresh(maxRes(), node->index, 1.0f);
		else
			refresh(maxRes(), node->index, 1.0f);

		// call recursively
		//if (endNode->neighbors.size() > 0)
			//drawNode(endNode);
	//}
}

//////////////////////////////////////////////////////////////////////
// set the intensity of the node
//////////////////////////////////////////////////////////////////////
void DAG::buildIntensity(NODE* root)
{
	if (root == NULL) return;

	// DAG segments
	int beginIndex = root->index;
	int begin[2];
	int end[2];
	float dWidth = _dx;
	float dHeight = _dy;

	for (unsigned int x = 0; x < root->neighbors.size(); x++)
	{
		static int count = 0;
		count++;
		if (count == 1) cout << "Now working on intensity:" << endl;
		if (!(count % 200)) cout << count<<" ";

		// get end node
		NODE* endNode = root->neighbors[x];
		int endIndex = endNode->index;

		// DAG segments
		begin[1] = beginIndex / _xRes;
		begin[0] = beginIndex - begin[1] * _xRes;
		end[1] = endIndex / _xRes;
		end[0] = endIndex - end[1] * _xRes;

		// set color
		if (endNode->leader)
			endNode->intensity = _leaderIntensity;
		else
		{
			// find max depth of current channel
			if (endNode->maxDepthNode == NULL)
				findDeepest(endNode, endNode->maxDepthNode);

			int maxDepth = endNode->maxDepthNode->depth;

			// calc standard deviation
			float stdDev = -(float)(maxDepth * maxDepth) / (float)(log(_secondaryIntensity) * 2.0f);

			// calc falloff
			float eTerm = -(float)(endNode->depth) * (float)(endNode->depth);
			eTerm /= (2.0f * stdDev);
			eTerm = exp(eTerm) * 0.5f;
			endNode->intensity = eTerm;
		}

		// call recursively
		if (endNode->neighbors.size() > 0)
			buildIntensity(endNode);
	}
}

//////////////////////////////////////////////////////////////////////
// DAG the tree offscreen
//////////////////////////////////////////////////////////////////////
float*& DAG::drawOffscreen(int scale)
{
	// allocate buffer
	_width = _xRes * scale;
	_height = _yRes * scale;
	_scale = scale;
	if (_offscreenBuffer) delete[] _offscreenBuffer;
	_offscreenBuffer = new float[_width * _height];

	// wipe the buffer
	for (int x = 0; x < _width * _height; x++)
		_offscreenBuffer[x] = 0.0f;

	// recursively DAG the tree
	drawOffscreenNode(_root);

	return _offscreenBuffer;
}

//////////////////////////////////////////////////////////////////////
// DAG the DAG segments
//////////////////////////////////////////////////////////////////////
void DAG::drawOffscreenNode(NODE* root)
{
	if (root == NULL) return;

	// DAG segments
	int beginIndex = root->index;
	int begin[2];
	int end[2];
	float dWidth = _dx;
	float dHeight = _dy;

	for (unsigned int x = 0; x < root->neighbors.size(); x++)
	{
		// get end node
		NODE* endNode = root->neighbors[x];
		int endIndex = endNode->index;

		// get endpoints 
		begin[0] = beginIndex % _xRes * _scale;
		begin[1] = beginIndex / _xRes * _scale;
		end[0] = endIndex % _xRes * _scale;
		end[1] = endIndex / _xRes * _scale;

		// make sure the one with the smaller x comes first
		if (end[0] < begin[0])
		{
			int temp[] = { end[0], end[1] };
			end[0] = begin[0];
			end[1] = begin[1];

			begin[0] = temp[0];
			begin[1] = temp[1];
		}

		// rasterize
		drawLine(begin, end, endNode->intensity);

		// call recursively
		if (endNode->neighbors.size() > 0)
			drawOffscreenNode(endNode);
	}
}

//////////////////////////////////////////////////////////////////////
// rasterize the line
// I assume all the lines are purely horizontal, vertical or diagonal
// to avoid using something like Bresenham
//////////////////////////////////////////////////////////////////////
void DAG::drawLine(int begin[], int end[], float intensity)
{
	int scaledX = _xRes * _scale;

	// if it is a horizontal line
	if (begin[1] == end[1])
	{
		for (int x = begin[0]; x < end[0]; x++)
		{
			int index = x + end[1] * scaledX;
			if (intensity > _offscreenBuffer[index])
				_offscreenBuffer[index] = intensity;
		}
		return;
	}

	// if it is a vertical line
	if (begin[0] == end[0])
	{
		int bottom = (begin[1] > end[1]) ? end[1] : begin[1];
		int top = (begin[1] > end[1]) ? begin[1] : end[1];

		for (int y = bottom; y < top; y++)
		{
			int index = begin[0] + y * scaledX;
			if (intensity > _offscreenBuffer[index])
				_offscreenBuffer[index] = intensity;
		}
		return;
	}

	// else it is diagonal
	int slope = (begin[1] < end[1]) ? 1 : -1;
	int interval = end[0] - begin[0];
	for (int x = 0; x <= interval; x++)
	{
		int index = begin[0] + x + (begin[1] + x * slope) * scaledX;
		if (intensity > _offscreenBuffer[index])
			_offscreenBuffer[index] = intensity;
	}
}

//////////////////////////////////////////////////////////////////////
// find deepest depth from current root
//////////////////////////////////////////////////////////////////////
void DAG::findDeepest(NODE* root, NODE*& deepest)
{
	deepest = root;

	for (unsigned int x = 0; x < root->neighbors.size(); x++)
	{
		NODE* child = root->neighbors[x];
		NODE* candidate = NULL;
		findDeepest(child, candidate);
		if (candidate->depth > deepest->depth)
			deepest = candidate;
	}
}



float* DAG::refresh(int Res,int index_, float inten) {
	static long int count = 0;
	count++;

	if (count == 1) {
		points = new float[Res * Res];
		for (int i = 0; i < Res; i++) {
			for (int j = 0; j < Res; j++) {
				points[i*Res + j] = 0;
			}
		}
	}


	points[index_] = inten;

	string filename=to_string(count);
	ofstream outfile(filename);

	for (int i = 0; i < Res; i++) {
		for (int j = 0; j < Res; j++) {
			outfile << points[i*Res + j]<<" ";
		}
		outfile << endl;
	}
	outfile.close();

	return points;
}


float* DAG::drawAll(NODE* _root,int Res, float* points_) {
	static bool first_in = true;

	if (first_in) {
		points_ = new float[Res * Res];
		for (int i = 0; i < Res; i++) {
			for (int j = 0; j < Res; j++) {
				points_[i*Res + j] = 0;
			}
		}
		first_in = false;
	}


	if (_root == NULL) {
		return NULL;
	}
	

	// draw segments
	int Index = _root->index;
	points_[Index] = _root->intensity;
	//int begin[2];
	//int end[2];
	//float dWidth = _dx;
	//float dHeight = _dy;

	for (unsigned int x = 0; x < _root->neighbors.size(); x++)
	{
		// get end node
		NODE* endNode = _root->neighbors[x];
		int endIndex = endNode->index;

		// draw segments
		//begin[1] = beginIndex / _xRes;
		//begin[0] = beginIndex - begin[1] * _xRes;
		//end[1] = endIndex / _xRes;
		//end[0] = endIndex - end[1] * _xRes;

		
		points_[endIndex]= endNode->intensity;

		if (endNode->neighbors.size() > 0)
			drawAll(endNode,Res,points_);
	}

	return points_;
}


void DAG::Alltofile(float* points_, int Res) {
	string filename = "final";
	ofstream outfile(filename);

	for (int i = 0; i < Res; i++) {
		for (int j = 0; j < Res; j++) {
			outfile << points_[i*Res + j] << " ";
		}
		outfile << endl;
	}
	outfile.close();

	

	delete[] points_;

	//cout << "Completed! You can close this program now." << endl;
}