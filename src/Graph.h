
#ifndef _GRAPH_H
#define _GRAPH_H

//
// NOTE: I do not split some functions into the .cpp file because forward declarations
//       make nasty work for linkers.
//
#include <iostream>
using std::istream;

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;

#include <fstream>
using std::ostream;

#include <iomanip>
using std::setw;

#include <list>
using std::list;

#include <vector>
using std::vector;

#include <climits>

#include <cassert>
#include <cmath>

// forward declarations so the compiler doesn't flip out
class Graph;
class Vertex;

const double UNINIT = 0.0;

class Edge
{
	friend class Graph;
public:
	Edge( Vertex* source, Vertex* destination, double weight = 0.0, double length = 1.0 )
	{
		assert(source); assert(destination); assert(weight>=0.0); assert(length>=1.0);
		src = source;
		dest = destination;
		this->weight = weight;
		this->length = length;
	}
	Vertex* getSource() { return src; }
	Vertex* getDestination() { return dest;	}
	double getWeight() const { return weight; }
	void setWeight(double Weight) { assert(Weight>=0.0); weight = Weight; }

	double getLength() const { return length; }
	void setLength(double Length) { assert(Length>=0); length = Length; }

	~Edge() { }
private:
	Vertex* src;
	Vertex* dest;
	double weight;
	double length;
};

class Vertex
{
	enum {IN_EDGE,OUT_EDGE};

	friend class Graph;
public:
	typedef list<Edge*>::const_iterator edge_const_iterator;
	typedef list<Edge*>::iterator edge_iterator;

	Vertex() { }

	Vertex( int Id, double Weight )
	{
		id = Id;
		weight = Weight;
	}

	int getID() const { return id; }

	~Vertex()
	{
		list<Edge*>::iterator itr =out_edgelist.begin(), end=out_edgelist.end();
		for ( ; itr != end; ++itr )
		{
			Edge* outEdgeToDelete = *itr;
			Vertex* neighbor = (*itr)->getDestination();

			neighbor->in_edgelist.remove( *itr );

			// look at this neighbor's out edge list, and see if we find our vertex.
			list<Edge*>::iterator neighborItr = neighbor->out_edgelist.begin(), neighborEnd=neighbor->out_edgelist.end();
			while ( neighborItr != neighborEnd && this != (*neighborItr)->getDestination() )
			{
				neighborItr++;
			}
			if ( neighborItr != neighborEnd )
			{

				Edge* toDelete = *neighborItr;
				neighbor->out_edgelist.erase( neighborItr );
				delete toDelete;
			}
			delete outEdgeToDelete;
		}
		out_edgelist.clear();
		in_edgelist.clear();
	}

	edge_iterator deleteInEdgeRef( edge_iterator& edge )
	{
		return deleteEdge( edge,IN_EDGE );
	}
	edge_iterator deleteOutEdgeRef( edge_iterator& edge )
	{
		return deleteEdge( edge,OUT_EDGE );
	}
	edge_iterator deleteEdge( edge_iterator& edge, int type )
	{
		edge_iterator next = edge;
		next++;

		Edge* toDelete = *edge;
		if ( type == IN_EDGE ) // deleting an IN edge
		{
			// remove it from the list
			//in_edgelist.erase(edge);

			// find it in the source and remove it from their out_edgelist
			Vertex* v = toDelete->getSource();
			for ( edge_iterator itr = v->outEdgesBegin(), end=v->outEdgesEnd();
					itr != end;
					itr++ )
			{
				if ( *itr == toDelete )
				{
					v->out_edgelist.erase(itr);
					break;
				}
			}
			delete toDelete;
		}
		else // deleting an OUT edge
		{
			//out_edgelist.erase(edge);


			Vertex* v = toDelete->getDestination();
			for ( edge_iterator itr = v->inEdgesBegin(), end=v->inEdgesEnd();
					itr != end;
					itr++ )
			{
				if ( *itr == toDelete )
				{
					v->in_edgelist.erase(itr);
					break;
				}
			}
			// delete the edge even if we do not find it in the other vertex
			delete toDelete;
		}
		return next;
	}




	double getWeight() const
	{
		return weight;
	}
	void setWeight( double weight )
	{
		this->weight = weight;
	}
	inline void insertOutEdge( Edge* e )
	{
		assert(e != NULL);
		out_edgelist.push_back( e );
	}
	inline void insertInEdge( Edge* e )
	{
		assert( e!=NULL );
		in_edgelist.push_back(e);
	}
	string toString() const
	{
		ostringstream oStr;
		oStr << "ID=" << id << "(" << weight << ")";
		return oStr.str();
	}

	edge_const_iterator outEdgesBegin() const { return out_edgelist.begin(); }
	edge_const_iterator outEdgesEnd() const { return out_edgelist.end(); }
	edge_iterator outEdgesBegin() { return out_edgelist.begin(); }
	edge_iterator outEdgesEnd() { return out_edgelist.end(); }

	edge_const_iterator inEdgesBegin() const { return in_edgelist.begin(); }
	edge_const_iterator inEdgesEnd() const { return in_edgelist.end(); }
	edge_iterator inEdgesBegin() { return in_edgelist.begin(); }
	edge_iterator inEdgesEnd() { return in_edgelist.end(); }
private:

	edge_iterator moveToHead( list<Edge*>& list, edge_iterator& itr )
	{
		list.push_front( *itr );
		return list.erase(itr);
	}

	edge_iterator moveOutEdgeToHead( edge_iterator& itr )
	{
		return moveToHead( out_edgelist, itr );
	}
	edge_iterator moveInEdgeToHead( edge_iterator& itr )
	{
		return moveToHead( in_edgelist, itr );
	}

	int id;
	double weight;

	list< Edge* > out_edgelist;
	list< Edge* > in_edgelist;
};




class Graph
{
public:
	typedef Vertex::edge_const_iterator edge_const_iterator;
	typedef Vertex::edge_iterator edge_iterator;


	Graph( int initDimension, double vertexInitialWeight, double edgeInitialWeight, double edgeInitialLength  );
	Graph( istream& readFromFile );
	~Graph();

	// addVertex( srcId, weight ) -- this will add a vertex to the graph with the identifier
	//      of ID. NOTE: if this vertex ID was already inserted, it will be overwritten.
	// PRE: 0 <= id < nVertices
	// POST: A new vertex with weight w is added into its slot within the vertices array
	void addVertex( int id, double weight = UNINIT );

	// addEdge( src, dest, weight ) -- this will add an directed edge (src,dest) of weight W
	//     to the graph. Note that the graph supports multiple edges between vertices.
	// PRE: src,dest on [0,nVertices) ; NO self-loops
	// POST: edge (src,dest) of weight W is inserted
	void addDirectedEdge( int srcId, int destId, double weight = UNINIT, double capacity = 1.0 );
	void addUndirectedEdge( int srcId, int destId, double weight = UNINIT, double capacity = 1.0 );

	// Vertex* getVertex(idx) -- get a pointer to a vertex with ID=idx.
	//    NOTE: it is up to the user of the class to determine if the vertex has already been inserted.
	// PRE: idx on [0,nVertices)
	// POST: the vertex ptr is returned
	Vertex* getVertex( int idx );

	// string toString() -- get a string representation of the graph
	// NOTE: the user of the class SHOULD insert all vertices before calling
	//    this. (Or junk will be returned.)
	string toString() const;
	string getDOTformat() const;

	static string getVertexDOT( int id, int inVertexTransition )
	{
		ostringstream out;
		out << "N" << id%inVertexTransition << ( (id<inVertexTransition) ? "in" : "out" );
		return out.str();
	}
	double getMinEnergy()
	{
		double min = INT_MAX;
		for ( int i = 0; i < nVertices; ++i )
		{
			double curr = getVertex(i)->getWeight();
			if ( curr < min )
			{
				min = curr;
			}
		}
		return min;
	}
	double getTotalEnergy()
	{
		double e = 0.0;
		for ( int i = 0; i < nVertices; ++i )
		{
			e += getVertex(i)->getWeight();
		}
		return e;
	}

	vector<double> getInNodeVertexWeight()
	{
		vector<double> inNode_energies( nVertices/2, 0 );
		if ( nVertices >= 1 )
		{
			for ( int i = 0; i < nVertices/2; ++i )
			{
				inNode_energies[i] = getVertex(i)->getWeight();
			}
		}
		return inNode_energies;
	}
	double dijkstra( int startId, int destinationId,
			list<Edge*>& pathOfEdges )
	{
		assert( startId >= 0 && startId < nVertices );
		assert( destinationId >= 0 && destinationId < nVertices );

		vector<double> vertexWeight( nVertices, INT_MAX );

		vector<Edge*> predEdge( nVertices, NULL );

		list<Vertex*> visited;

		vertexWeight[startId] = 0.0;

		visited.push_front( &vertices[startId] );

		// while we haven't found the destination and there are unchecked nodes, continue
		while ( !visited.empty() )
		{
			Vertex* smallest = visited.front();
			visited.pop_front();

			// if we have have found the destination, STOP.
			if ( smallest->getID() == destinationId )
			{
				visited.clear();
				break;
			}

			// visit all adjacent nodes of the vertex with the smallest weight
			edge_iterator edgeItr = smallest->outEdgesBegin(), end = smallest->outEdgesEnd();
			for( ; edgeItr != end; ++edgeItr )
			{
				Edge* theEdge = *edgeItr;

				// vertexWeight[smallest->get()] is the length from the SOURCE node to the 'smallest' node.
				// IF we can cross the edge (hasn't used its E uses) and
				//   we can visit the dest node (it isn't dead) and
				//     the length from SOURCE to Smallest->neighbor, passing through Smallest, is less than the
				//     previous length then save the new length
				// MUST BE _STRICTLY_ LESS THAN
				if ( theEdge->getWeight() >= 1.0 &&
						theEdge->dest->getWeight() >= 1.0 &&
						( vertexWeight[smallest->getID()] + theEdge->getLength() < vertexWeight[theEdge->dest->getID()] ))
				{
					// if we have NOT seen the vertex before then we need to add it to the set
					if ( vertexWeight[theEdge->dest->getID()] == INT_MAX )
					{
						visited.push_back( theEdge->dest );
					}

					vertexWeight[theEdge->dest->getID()] = vertexWeight[smallest->getID()] + theEdge->getLength();
					predEdge[theEdge->dest->id] = theEdge;
				}
			}

			// find smallest and place at front
			if ( !visited.empty() )
			{
				list<Vertex*>::iterator itr = visited.begin(), end=visited.end(), loc = itr;
				for ( ; itr != end; ++itr )
				{
					if ( vertexWeight[(*itr)->getID()] < vertexWeight[(*loc)->getID()] )
						loc = itr;
				}
				smallest = *loc;
				visited.erase( loc );
				visited.push_front( smallest );

				assert( smallest->getWeight() >= 1.0 );
			}
		}

		int theId = destinationId;


		if ( predEdge[destinationId] == NULL )
		{
			return -1;
		}


		while ( theId != startId )
		{
			assert( predEdge[theId]->getLength() >= 1.0 );
			assert( predEdge[theId]->dest->getWeight() >= 1.0 );

			pathOfEdges.push_front( predEdge[theId] );
			theId = predEdge[theId]->src->getID();
		}

		return vertexWeight[destinationId];
	}

	string toStrinSimple() const
	{
		// NodeID NodeWeight <EdgeDest EdgeWeight EdgeLength>* \n
		// int    double     <int      double     double>*
		ostringstream out;
		out << nVertices << "\n";
		for ( int v = 0; v < nVertices; ++v )
		{
			out <<  vertices[v].id << " " << vertices[v].weight;
			edge_const_iterator edgeItr = vertices[v].outEdgesBegin(),
				end = vertices[v].outEdgesEnd();

			for ( ; edgeItr != end; ++edgeItr )
			{
				out << " " << (*edgeItr)->dest->id
					<< " " << (*edgeItr)->weight
					<< " " << (*edgeItr)->length;
			}
			out << "\n";
		}

		return out.str();
	}

	string getAllNodeWeights() const
	{
		ostringstream out;
		int dim = (int)sqrt( (double)nVertices );
		out << dim << '\n';
		for ( int i = 0; i < dim; ++i )
		{
			for ( int j = 0; j < dim; ++j )
			{
				out << vertices[i*dim+j].getWeight() << ' ';
			}
			out << '\n';
		}
		return out.str();
	}

protected:
	int nVertices;
	Vertex* vertices;
};

#endif
