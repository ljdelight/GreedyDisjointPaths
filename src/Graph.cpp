#include "Graph.h"

#include <vector>
using std::vector;

#include <iostream>
using std::cerr;
using std::istream;
using std::cout;

#include <sstream>
using std::ostringstream;

#include <new>
using std::nothrow;

#include <iomanip>
using std::setw;
using std::setprecision;
using std::fixed;

#include <cstdlib>
#include <cmath>
#include <climits>


Graph::Graph( int initDimension, double vertexInitialWeight, double edgeInitialWeight, double edgeInitialLength )
{
	assert( initDimension > 0);

	nVertices = 2* initDimension * initDimension;
	int nVert = nVertices / 2;

	vertices = new (nothrow) Vertex[nVertices];

	// initialize all of the graph's vertices
	for ( int i = 0; i < nVertices; ++i )
	{
		addVertex(i,vertexInitialWeight );
	}


	for (int i = 0; i < initDimension; ++i)
	{
		for (int j = 0; j < initDimension; ++j)
		{
			int id = i * initDimension + j;

			// THE 'IN'  NODE IS id
			// THE 'OUT' NODE IS id+nVert

			// connect 'IN' to 'OUT'
			addDirectedEdge( id,id+nVert, edgeInitialWeight, edgeInitialLength );


			// add an edge to the left, if we're not on the first column
			if (j != 0)
			{
				addDirectedEdge(id + nVert, id - 1, edgeInitialWeight,
						edgeInitialLength);
			}
			// add an edge to the right, if we're not on the last column
			if ( j != initDimension-1 )
			{
				addDirectedEdge(id + nVert, id + 1, edgeInitialWeight,
									edgeInitialLength);
			}
			// add an edge UP, if we're not on the first row
			if (i != 0)
			{
				addDirectedEdge(id + nVert , (i - 1) * initDimension + j,
						edgeInitialWeight, edgeInitialLength);
			}
			// add an edge DOWN if we're not on the last row
			if ( i != initDimension-1 )
			{
				addDirectedEdge(id + nVert , (i + 1) * initDimension + j,
									edgeInitialWeight, edgeInitialLength);
			}
		}
	}
	if ( vertices == NULL )
	{
		cerr << "ERROR: could not allocate space for vertices in Graph::Graph(int). Exiting.\n";
		exit(1);
	}
}

Graph::Graph( istream& inputFileStream )
{
	if ( !inputFileStream )
	{
		cerr << "ERROR: Graph::Graph istream& maked as not 'good'.\n";
		throw;
	}
	int dim = 0;

	if ( ! (inputFileStream>>dim) )
	{
		cerr << "ERROR: Graph::Graph could not read the dimension from file!\n";
		throw;
	}
	nVertices = dim*dim;

	vertices = new (nothrow) Vertex[nVertices];
	if ( vertices == NULL )
	{
		cerr << "ERROR: could not allocate space for vertices in Graph::Graph(int). Exiting.\n";
		exit(1);
	}

	double nodeWeight;
	for ( int i = 0; i < dim; ++i )
	{
		for ( int j = 0; j < dim; ++ j )
		{
			int id = i*dim+j;

			if ( inputFileStream >> nodeWeight )
			{
				addVertex( id, nodeWeight );

				// add an edge to the left, if we're not on the first column
				if ( j != 0 )
				{
					addUndirectedEdge( id, id -1, nodeWeight );
				}
				// add an edge UP, if we're not on the first row
				if ( i != 0 )
				{
					addUndirectedEdge( id, (i-1)*dim+j, nodeWeight );
				}
			}
			else
			{
				cerr << "ERROR: Graph::Graph did not read a " << dim << " by " << dim << " network. Validate input.\n";
				throw;
			}
		}
	}
	if ( (inputFileStream >> nodeWeight) )
	{
		cerr << "WARNING: Graph::Graph finished reading the graph but there remain unread node weights on the input stream.\n";
		throw;
	}

}

Graph::~Graph()
{
	// this will call each respective Edge dtor, freeing all edges.
	for ( int i = 0; i < nVertices; ++i )
	{
		vertices[i].in_edgelist.clear();
	}
	delete [] vertices;
	vertices = NULL;
}


// addVertex( srcId, weight ) -- this will add a vertex to the graph with the identifier
//      of ID. NOTE: if this vertex ID was already inserted, it will be overwritten.
// PRE: 0 <= id < nVertices
// POST: A new vertex with weight w is added into its slot within the vertices array
void Graph::addVertex( int id, double weight )
{
	assert( id >= 0 ); assert( id < nVertices );

	vertices[id] = Vertex( id, weight );
}

// addEdge( src, dest, weight ) -- this will add an directed edge (src,dest) of weight W
//     to the graph. Note that the graph supports multiple edges between vertices.
// PRE: src,dest on [0,nVertices) ; NO self-loops
// POST: edge (src,dest) of weight W is inserted
void Graph::addDirectedEdge( int srcId, int destId, double weight, double length )
{
	// the vertices must EXIST! Also no loops.
	assert( srcId != destId );
	assert( srcId >= 0 &&  srcId < nVertices );
	assert( destId >= 0 && destId < nVertices );

	Vertex *src  = &vertices[srcId],
		   *dest = &vertices[destId];

	// create the edge and add it to the source vertex.
	Edge* srcToDest = new (nothrow) Edge( src, dest, weight, length );

	if ( srcToDest == NULL )
	{
		cerr << "ERROR: Could not allocate Edge in Graph::addEdge(int,int,double). Exiting.\n";
		exit(1);
	}

	// the sending and receiving vertex inserts this edge into their in, out edge list respectively.
	src->insertOutEdge(srcToDest);
	dest->insertInEdge(srcToDest);
}

void Graph::addUndirectedEdge( int srcId, int destId, double weight, double length )
{
	addDirectedEdge(srcId,destId,weight,length);
	addDirectedEdge(destId,srcId,weight,length);
}


string Graph::getDOTformat() const
{
	ostringstream out;
	int inVertexTransition = nVertices/2;
	for ( int v = 0; v < nVertices; ++v )
	{
		int id = vertices[v].id;
		if ( v < inVertexTransition )
		{

			cout << "\t\"N" << id << "in\""  << " [ label=< N<sub>" << id << "</sub><sup>in</sup> "
				 << vertices[v].weight << ">, color=red ];\n";
			cout << "\t\"N" << id << "out\"" << " [ label=< N<sub>" << id << "</sub><sup>out</sup> "
				 << vertices[v+inVertexTransition].weight << ">, color=blue ];\n";
		}

		edge_const_iterator edgeItr = vertices[v].outEdgesBegin(),
			end = vertices[v].outEdgesEnd();

		for ( ; edgeItr != end; ++edgeItr )
		{
			int destId = (*edgeItr)->getDestination()->getID();

			out << '"' << getVertexDOT(id,inVertexTransition) << '"'
				<< " -> "
				<< '"' << getVertexDOT(destId,inVertexTransition) << '"'
				<< " [ label=<" << (*edgeItr)->weight << "> ]\n";
		}
	}

	return out.str();
}


string Graph::toString() const
{
	ostringstream out;

	out << ":::ID=? (nodeWeight=?) Neighbors: ID(edgeWeightToNeighbor,length)\n";
	for ( int v = 0; v < nVertices; ++v )
	{
		out << "ID=" << vertices[v].id << " (W=" << vertices[v].weight << ")";
		edge_const_iterator edgeItr = vertices[v].outEdgesBegin(),
			end = vertices[v].outEdgesEnd();

		out << " OUT: {";
		// print the out edges
		if ( edgeItr != end )
		{
			out << (*edgeItr)->dest->id;
			out << "(" << (*edgeItr)->weight << "," << (*edgeItr)->length << ")";
			edgeItr++;
		}
		for ( ; edgeItr != end; ++edgeItr )
		{
			out << "," << (*edgeItr)->dest->id;
			out << "(" << (*edgeItr)->weight << "," <<  (*edgeItr)->length << ")";
		}


		edgeItr = vertices[v].inEdgesBegin();
		end = vertices[v].inEdgesEnd();
		out << "}  IN: {";

		// print the IN edges
		if ( edgeItr != end )
		{
			out << (*edgeItr)->src->id;
			out << "(" << (*edgeItr)->weight << "," << (*edgeItr)->length << ")";
			edgeItr++;
		}
		for ( ; edgeItr != end; ++edgeItr )
		{
			out << "," << (*edgeItr)->src->id;
			out << "(" << (*edgeItr)->weight << "," << (*edgeItr)->length << ")";
		}
		out << "}\n";
	}

	return out.str();
}


// Vertex* getVertex(idx) -- get a pointer to a vertex with ID=idx.
//    NOTE: it is up to the user of the class to determine if the vertex has already been inserted.
// PRE: idx on [0,nVertices)
// POST: the vertex ptr is returned
Vertex* Graph::getVertex( int id )
{
	assert( id >= 0 ); assert( id < nVertices );
	return &vertices[id];
}
