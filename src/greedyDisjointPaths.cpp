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
			<< "    OR stdin=[<dim> <E> < <src> <dest> >+]" << endl;
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
	int dimension = 0, m = 0, readRequests=0, satisfiedRequests = 0;
	double E = 0.0, beta = 0.0, res = -1.0;

	if ( argc==2  )
	{
		out.close();
		writeToFile = true;
		out.open(argv[1]);
		if ( !out  )
		{
			cerr << "ERROR: could not open " << argv[1] << endl;
			exit(1);
		}
	}
	if ( argc==5  )
	{
		out.close();
		writeToFile = true;
		out.open(argv[4]);
		if ( !out  )
		{
			cerr << "ERROR: could not open " << argv[4] << endl;
			exit(1);
		}
		argc -= 1;
	}
	if ( argc==4 )
	{
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
	}
	else
	{
		// read dim, initial network energy, and total number of messages
		IN() >> dimension >> E;
		if ( !IN() ) { cerr << "ERROR: cannot read dim, energy; invalid input format." << endl; exit(1); }

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

			//OUT() << g.toString()<< endl;
		int x,y;
		while ( IN_PAIR() >> x >> y )
		{
			readRequests += 1;
			pathOfEdges.clear();
			res = g.dijkstra( x, y, pathOfEdges );

			if ( res != -1 )
			{
				satisfiedRequests += 1;

				traversePathAndUpdate(beta, pathOfEdges);


				//cout << "\n}\n";

				//exit(0);

			}
			else
			{
				//OUT() << "CANNOT ROUTE " << x << " TO " << y << endl;
			}
			OUT() << g.toStrinSimple();
			OUT() << x << " " << y;
			list<Edge*>::const_iterator itr = pathOfEdges.begin(),
					end = pathOfEdges.end();
			for ( ; itr != end; ++itr )
			{
				OUT() << " " << (*itr)->getSource()->getID() << " " << (*itr)->getDestination()->getID();
			}

			OUT() << endl << "@";
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
