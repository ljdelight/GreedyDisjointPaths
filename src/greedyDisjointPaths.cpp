// Lucas Burson
//
// Implementation of GDP algorithm in Section V part B.
//
// TODO: Allow for param to disable output that is destined for the Graphviz builder.
//
// NOTE: Here I show how we compute 'm' using only the dimension. Much easier than
//       tracking |E| and |V| for graph G or G'. Assume square graphs.
// m = |E(G')| = 2*|E(G)|+|V(G)|
//             = 2(2*dim*(dim-1)) + dim^2
//             = 4dim^2 - 4dim + dim^2
//             = 5dim^2 - 4 dim
// E = initial energy of all nodes
// \beta = m^{1/(E+1)} > 1
//
#include "Graph.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;


bool readFromFile=false;
ifstream in;
bool writeToFile=false;
ofstream out("graphviz_intermed_file.mine");
ostream& OUT()
{
	return out;
}
istream& IN()
{
	return cin;
}
istream& IN_PAIR()
{
	if ( !readFromFile )
	{
		return cin;
	}
	else
	{
		return in;
	}
}

void USAGE()
{
	cerr << "programName [<Dim> <Energy_FILE_NAME> <Msg_Pairs_FILE_NAME>]\n"
			<< "  Energy_file contains an int representing the initial energy of the entire network\n"
			<< "  Msg_Pairs contains pairs of ints separated by whitespace; first int routes to second int"<< endl;
	exit(1);
}

void printEdgeListInfo( const list<Edge*>& pathOfEdges, ostream& outStream );
void printPath( const list<Edge*>& pathOfEdges, ostream& outStream );
void traversePathAndUpdate( double beta, list<Edge*>& pathOfEdges );


int main( int argc, char* argv[] )
{
	srand( (unsigned int)time(NULL) );

	if ( argc!=4 ) 
	{
		USAGE();
	}

	//
	// m = |E(G')| = 2*|E(G)|+|V(G)|
	//             = 2(2*dim*(dim-1)) + dim^2
	//             = 4dim^2 - 4dim + dim^2
	//             = 5dim^2 - 4 dim
	// E = initial energy of all nodes
	// \beta = m^{1/(E+1)} > 1
	//

	list<Edge*> pathOfEdges;
	int dimension = 0, m = 0, satisfiedRequests = 0;
	double E = 0.0, beta = 0.0, res = -1.0;

	readFromFile = true;
	dimension = atoi( argv[1] );
	in.open(argv[2]);
	if ( !in  )
	{
		cerr << "ERROR: could not open " << argv[2] << endl;
		exit(1);
	}
	in >> E;
	in.close();

	in.open(argv[3]);
	if ( !in )
	{
		cerr << "ERROR: could not open " << argv[3] << endl;
		exit(1);
	}



	if ( !(dimension>=2) ) { cerr << "ERROR: dimension must be >=2" << endl; exit(1); }
	if ( !(E>=1) ) { cerr << "ERROR: initial energy must be >=1" << endl; exit(1); }

	m = 5*dimension*dimension - 4*dimension;
	beta = pow(m, 1.0/(E+1) );

	try
	{
		Graph g( dimension, E, E, 1.0 );

		OUT() << g.toStrinSimple() << "@\n";

		double initEnergy = g.getTotalEnergy();
		int x,y;
		list< pair<int,int> > unsatisfiedRequests;

		// read all the possible requests
		while ( IN_PAIR() >> x >> y )
		{
			unsatisfiedRequests.push_back( make_pair(x,y) );
		}
		OUT() << "attempting to route " << unsatisfiedRequests.size() << " requests" << endl;
		assert( !unsatisfiedRequests.empty() );




		pair<int,int> min_pair;
		list<Edge*> min_pathOfEdges;
		double min_shortestPath = INT_MAX;


		// while we have unsatisfied requests,
		while ( !unsatisfiedRequests.empty() )
		{
			// initialize these
			min_pair = make_pair(-1,-1);
			min_pathOfEdges.clear();
			min_shortestPath = INT_MAX;


			// For each unsatisfied request, find its shortest-path.
			// Select the minimum-weight shortest path. Remove this path from the unsatisfiedRequest list.
			// Traverse this path
			list< pair<int,int> >::iterator itr = unsatisfiedRequests.begin(),
								end = unsatisfiedRequests.end();
			while ( itr != end )
			{
				pathOfEdges.clear();
				res = g.dijkstra( (*itr).first, (*itr).second, pathOfEdges );

				// remove points that we can't satisfy
				if ( res == -1 )
				{
					OUT() << "could not connect the request (" << (*itr).first << ", " << (*itr).second <<")" << endl;
					itr = unsatisfiedRequests.erase( itr );
				}
				// if we find a better shortest-path, save it.
				else if ( res < min_shortestPath )
				{
					OUT() << "Round: Better shortest path: (" << (*itr).first << "," << (*itr).second << ") of length " << res << endl;
					min_pair = *itr;
					min_pathOfEdges = pathOfEdges;
					min_shortestPath = res;

					itr++;
				}
				else
				{
					itr++;
				}
			}


			if ( min_shortestPath != INT_MAX )
			{
				assert( min_pair.first != -1 );assert( min_pair.second != -1 );

				OUT () <<"Satisfied (" << min_pair.first << "," << min_pair.second << ")" << endl;
				unsatisfiedRequests.remove( min_pair );
				satisfiedRequests += 1;
				traversePathAndUpdate(beta, min_pathOfEdges);
			}
		}


		double finalEnergy = g.getTotalEnergy();

		cout << g.getMinEnergy() << " ";
		cout << (initEnergy - finalEnergy)/2.0 << " ";
		cout << satisfiedRequests << " ";

		vector<double> inNode_energies = g.getInNodeVertexWeight();
		if ( inNode_energies.size() )
		{
			cout << inNode_energies[0];
			for ( unsigned int i = 1; i < inNode_energies.size(); ++i )
			{
				cout << "," << inNode_energies[i];
			}
		}

		cout << endl;


	} catch( ... )
	{
		cerr << "ERROR: An exception was thrown. Exiting." << endl;
		exit(1);
	}
	exit(0);
}




void printEdgeListInfo( const list<Edge*>& pathOfEdges, ostream& outStream )
{
	outStream << "(src->dest,weight,length)\n";
	for ( list<Edge*>::const_iterator itr = pathOfEdges.begin(), end = pathOfEdges.end();
				itr != end; ++itr )
	{
		Edge* p = *itr;
		outStream << "(" << p->getSource()->getID() << "->"
				         << p->getDestination()->getID() << ","
				         << p->getWeight() << ","
				         << p->getLength() << ") ";
	}
	outStream << "\n";
}
void printPath( const list<Edge*>& pathOfEdges, ostream& outStream )
{
	Edge* p = NULL;
	for ( list<Edge*>::const_iterator itr = pathOfEdges.begin(), end = pathOfEdges.end();
					itr != end; ++itr )
	{
		p = *itr;
		outStream << p->getSource()->getID() << "->";
	}
	if ( p != NULL )
	{
		outStream << p->getDestination()->getID() << "\n";
	}
}
void traversePathAndUpdate( double beta, list<Edge*>& pathOfEdges )
{
	Edge *p = NULL;
	Vertex *src = NULL;
	for ( list<Edge*>::const_iterator itr = pathOfEdges.begin(), end = pathOfEdges.end();
						itr != end; ++itr )
	{
		p = *itr;
		src = p->getSource();

		// length must always be positive
		assert( p->getLength() >= 1.0 );

		// weight can be 0, BUT we are traversing the edge so it must be crossable.
		assert( p->getWeight() >= 1.0 );

		src->setWeight( src->getWeight() - 1.0 );
		p->setLength( p->getLength() * beta );
		p->setWeight( p->getWeight() - 1.0 );
	}
}
