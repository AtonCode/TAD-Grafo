// =========================================================================
// @author: David Enrique Palacios Garc�a
// =========================================================================
#ifndef __Graph__hxx__
#define __Graph__hxx__

#include <cstdlib>
#include <queue>
#include <stack>
#include <limits>
#include <sstream>
#include <algorithm>    // std::reverse

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
Graph< _TValue, _TCost >::
Graph( )
{
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
Graph< _TValue, _TCost >::
~Graph( )
{
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
unsigned long Graph< _TValue, _TCost >::
AddVertex( const TValue& v )
{
    this->m_Vertices.push_back( v );
    return( this->m_Vertices.size( ) - 1 );
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
void Graph< _TValue, _TCost >::
SetArc( unsigned long i, unsigned long j, const TCost& c )
{
    this->m_Costs[ i ][ j ] = c;
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
unsigned long Graph< _TValue, _TCost >::
GetNumberOfVertices( ) const
{
    return( this->m_Vertices.size( ) );
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
const typename Graph< _TValue, _TCost >::
TValue& Graph< _TValue, _TCost >::
GetVertex( unsigned long i ) const
{
    return( this->m_Vertices[ i ] );
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
bool Graph< _TValue, _TCost >::
HasArc( unsigned long i, unsigned long j ) const
{
    typename TMatrix::const_iterator rIt = this->m_Costs.find( i );
    if( rIt != this->m_Costs.end( ) )
        return( rIt->second.find( j ) != rIt->second.end( ) );
    return( false );
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
const typename Graph< _TValue, _TCost >::
TCost& Graph< _TValue, _TCost >::
GetCost( unsigned long i, unsigned long j ) const
{
    static const _TCost inf_cost = std::numeric_limits< _TCost >::max( );
    if( this->HasArc( i, j ) )
        return( this->m_Costs[ i ][ j ] );
    else
        return( inf_cost );
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
void Graph< _TValue, _TCost >::
PrintPlaneGraph( ) const
{
    typename TVertices::const_iterator vIt = this->m_Vertices.begin( );
    for( ; vIt != this->m_Vertices.end( ); vIt++ )
        std::cout << *vIt << " ";
    std::cout << std::endl;
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
void Graph< _TValue, _TCost >::
PrintPreorderGraph( unsigned long i, std::vector<bool>& m )
{
    if(m[i]==false)
    {
        std::cout<<this->m_Vertices[i]<<" ";
        m[i]=true;
        int next;
        std::vector<unsigned long> adyacencia = AdyacenciasVertex(i);
        for(int k=0; k<adyacencia.size(); k++)
        {
            next = adyacencia[k];
            if(m[next]==false)
            {
                PrintPreorderGraph(next, m);
            }
        }
        adyacencia.clear();
    }
}

template< class _TValue, class _TCost >
std::vector< unsigned long > Graph< _TValue, _TCost >::
AdyacenciasVertex(unsigned long i)
{
    std::vector< unsigned long > adyacencias;
    typename TMatrix::const_iterator cIt = this->m_Costs.find( i );
    if( cIt != this->m_Costs.end( ) )
    {
        typename TRow::const_iterator rIt = cIt->second.begin( );
        for( ; rIt != cIt->second.end( ); rIt++ )
            adyacencias.push_back( rIt->first );
    }
    return adyacencias;
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
void Graph< _TValue, _TCost >::
PrintLevelsGraph( unsigned long i ) const
{
    std::vector< bool > m( this->m_Vertices.size( ), false );
    std::queue< unsigned long > q;
    q.push( i );
    while( !q.empty( ) )
    {
        unsigned long n = q.front( );
        q.pop( );

        if( m[ n ] ) //si ya tiene true
            continue;

        std::cout << this->m_Vertices[ n ] << " ";
        m[ n ] = true;

        typename TMatrix::const_iterator cIt = this->m_Costs.find( n );
        if( cIt != this->m_Costs.end( ) )
        {
            typename TRow::const_iterator rIt = cIt->second.begin( );
            for( ; rIt != cIt->second.end( ); rIt++ )
                q.push( rIt->first );

        } // fi

    } // elihw
    std::cout << std::endl;
}

// -------------------------------------------------------------------------
template< class _TValue, class _TCost >
void Graph< _TValue, _TCost >::
PrintGraphAsPNG( const std::string& filename ) const
{
    std::stringstream str;
    str << "echo \"digraph G{";

    typename TMatrix::const_iterator cIt = this->m_Costs.begin( );
    for( ; cIt != this->m_Costs.end( ); cIt++ )
    {
        typename TRow::const_iterator rIt = cIt->second.begin( );
        for( ; rIt != cIt->second.end( ); rIt++ )
        {
            str << cIt->first << "->" << rIt->first << " ";
            str
                    << cIt->first << " [label=\""
                    << this->m_Vertices[ cIt->first ] << "\"]; ";
            str
                    << rIt->first << " [label=\""
                    << this->m_Vertices[ rIt->first ] << "\"]; ";

        } // rof

    } // rof

    str << "}\" | dot -Tpng > " << filename;
    std::system( str.str( ).c_str( ) );
}
// -------------------------------------------------------------------------------------

template< class _TValue, class _TCost >
std::vector< unsigned long > Graph< _TValue, _TCost >::
Dijkstra( unsigned long i, unsigned long j ) const
{
    struct _TNode
    {
        unsigned long Vertex;
        unsigned long Parent;
        TCost         Cost;

        _TNode(unsigned long v, unsigned long p, TCost c) //Creador del nodo
        {
            this->Vertex = v;
            this->Parent = p;
            this->Cost = c;
        }

        bool operator<( const _TNode& right)
        {
            return ( right.Cost < this->Cost );
        }
    };


    //Inicializaci�n
    std::vector<bool> marks(this->m_Vertices.size(), false); //Vector de marcas de los nodos y se inicializan en false
    std::vector<long> t(this->m_Vertices.size(), -1); //Vector de los vertices y se inicializan en -1

    std::vector<_TNode> q; //Se crea una cola de prioridad
    q.push_back(_TNode(i,j,0)); //Agrega el nodo a la cola de prioridad

    std::make_heap(q.begin(), q.end()); //Organiza los puntos de la cola de prioridad
    while(!q.empty()) //Mientras que la cola de prioridad NO este vacia
    {
        std::pop_heap(q.begin(), q.end()); //Send min to back
        _TNode n = q.back(); //Get min
        q.pop_back(); //Borra el nodo de la cola de prioridad

        if ( marks[ n.Vertex ] ) //Si el nodo ya est� marcado, sale y continua
        {
            continue;
        }

        marks[n.Vertex] = true; //Marca la visita del nodo
        t[ n.Vertex ]=n.Parent;//A�ade un padre


        typename TMatrix::const_iterator cIt = this->m_Costs.find( n.Vertex);
        if(cIt != this->m_Costs.end( ))
        {
            typename TRow::const_iterator rIt = cIt->second.begin( );
            for( ; rIt != cIt->second.end( ); rIt++)
            {
                q.push_back( _TNode(rIt->first, n.Vertex, rIt->second + n.Cost));
                std::push_heap(q.begin(), q.end( ));
            } //end for
        } //end if

    } //end while

    //Backtracking
    std::vector<unsigned long> camino;
    if(marks[j])
    {
        unsigned long k=j; //k es la llegada
        while(i!=k) //mientras el hijo y el padre no sean iguales (i es raiz)
        {
            camino.push_back(k);
            k=t[k];
        }
        camino.push_back(k);
    }
    std::reverse(camino.begin(), camino.end());
    return camino;
}

//------------------------------------------------------------------------------------------------
template< class _TValue, class _TCost >
std::vector < std::vector < unsigned long > > Graph< _TValue, _TCost >::
Mapaavector()
{
    std::vector < std::vector < unsigned long > > res;
    for(int i=0; i<this->m_Vertices.size(); i++)
    {
        std::vector<unsigned long> fila;
        for(int j=0; j<this->m_Vertices.size(); j++)
        {
            fila.push_back(this->m_Costs[i][j]);
        }
        res.push_back(fila);
    }
    return res;
}

//------------------------------------------------------------------------------------------------
template <class _TValue, class _TCost>
std::vector < std::vector < unsigned long > > Graph<_TValue, _TCost>::
FloydWarshall()
{
    int cn = this->GetNumberOfVertices(); //cantidad de nodos
// Devuelve una matriz con las distancias m�nimas de cada nodo al resto de los v�rtices.

    std::vector< std::vector<unsigned long> > path = this->Mapaavector(); //convierte el mapa a un vector

    for(int i = 0; i < cn; i++)
        path[i][i] = 0; //diagonales en 0


    for(int i=0 ; i<this->GetNumberOfVertices(); i++)
    {
        for(int j=0; j<this->GetNumberOfVertices(); j++ )
            if(path[i][j]==0)
                path[i][j]=99999;
    }

    for(int k = 0; k < cn; k++)
        for(int i = 0; i < cn; i++)
            for(int j = 0; j < cn; j++)
            {
                int dt = path[i][k] + path[k][j];
                if(path[i][j] > dt)
                    path[i][j] = dt;
            }
    return path;
}

//------------------------------------------------------------------------------------------------

template <class _TValue, class _TCost>
std::deque< _TValue > Graph<_TValue, _TCost>::
getFloyd( std::vector < std::vector < unsigned long > > iaPath, unsigned long iSrc, unsigned long iDest)
{
    std::deque< _TValue > camino;
    std::stack<unsigned long> S;
    unsigned long iFrom;

    S.push(iDest);

    iFrom = iaPath[iSrc][iDest];//Get parent or from vertex ID.
    while(1)
    {
        if(iFrom == 99999 ) //Check path exists.
        {
            std::cout << "There is no path\n";
            return camino;
        }
        S.push(iFrom);

        if( iSrc == iFrom )
            break; //Check if from vertex and source vertex are the same. Is so, touchdown!!!
        iFrom = iaPath[iSrc][iFrom];
    }

//cout << iSrc << "==>" << endl;
    while(!S.empty())
    {
        std::cout << S.top() << "==>";
        camino.push_back(this->GetVertex(S.top()));
        S.pop();
    }
    std::cout << std::endl;


    return camino;
}

//------------------------------------------------------------------------------------------------
template <class _TValue, class _TCost>
bool Graph<_TValue, _TCost>::
isSafe(unsigned long v, unsigned long path[], unsigned long pos)
{
    /* Check if this vertex is an adjacent
    vertex of the previously added vertex. */
    if (this->m_Costs [path[pos - 1]][ v ] == 0)
        return false;

    /* Check if the vertex has already been included.
    This step can be optimized by creating
    an array of size V */
    for (int i = 0; i < pos; i++)
        if (path[i] == v)
            return false;

    return true;
}
//------------------------------------------------------------------------------------------------
template <class _TValue, class _TCost>
bool Graph<_TValue, _TCost>::
hamCycleUtil(unsigned long path[], unsigned long pos)
{
    unsigned long V = this->GetNumberOfVertices();
    /* base case: If all vertices are
    included in Hamiltonian Cycle */
    if (pos == V)
    {
        // And if there is an edge from the
        // last included vertex to the first vertex
        if (this->m_Costs[path[pos - 1]][path[0]] == 1)
            return true;
        else
            return false;
    }

    // Try different vertices as a next candidate
    // in Hamiltonian Cycle. We don't try for 0 as
    // we included 0 as starting point in hamCycle()
    for (int v = 1; v < V; v++)
    {
        /* Check if this vertex can be added
        // to Hamiltonian Cycle */
        if (isSafe(v, path, pos))
        {
            path[pos] = v;

            /* recur to construct rest of the path */
            if (hamCycleUtil (path, pos + 1) == true)
                return true;

            /* If adding vertex v doesn't lead to a solution,
            then remove it */
            path[pos] = -1;
        }
    }

    /* If no vertex can be added to
    Hamiltonian Cycle constructed so far,
    then return false */
    return false;
}

//------------------------------------------------------------------------------------------------
template <class _TValue, class _TCost>
std::deque <unsigned long> Graph<_TValue, _TCost>::
hamCycle()
{
    std::deque < unsigned long > res;
    unsigned long V = this->GetNumberOfVertices();
    unsigned long *path = new unsigned long[V];
    for (int i = 0; i < V; i++)
        path[i] = -1;

    /* Let us put vertex 0 as the first vertex in the path.
    If there is a Hamiltonian Cycle, then the path can be
    started from any point of the cycle as the graph is undirected */
    path[0] = 0;

    if (hamCycleUtil(path, 1) == true )
    {
       for(unsigned long i=0; i<V; i++)
        res.push_back(path[i]);
       res.push_back(path[0]);
    }
    return res;
}






#endif // __Graph__hxx__

// eof - Graph.hxx
