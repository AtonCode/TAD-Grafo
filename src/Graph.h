// =========================================================================
// @author: David Enrique Palacios Garc�a
// =========================================================================

#ifndef __Graph__h__
#define __Graph__h__

#include <map>
#include <string>
#include <vector>
#include <deque>

template< class _TValue, class _TCost >
class Graph
{
public:
  typedef _TValue TValue;
  typedef _TCost  TCost;

  typedef std::vector< TValue >            TVertices;
  typedef std::map< unsigned long, TCost > TRow;
  typedef std::map< unsigned long, TRow >  TMatrix;

public:
  Graph( );
  virtual ~Graph( );

  unsigned long AddVertex( const TValue& v );
  void SetArc( unsigned long i, unsigned long j, const TCost& c );
  std::vector< unsigned long> Dijkstra( unsigned long i, unsigned long j ) const;

  unsigned long GetNumberOfVertices( ) const;
  const TValue& GetVertex( unsigned long i ) const;

  bool HasArc( unsigned long i, unsigned long j ) const;
  const TCost& GetCost( unsigned long i, unsigned long j ) const;

  std::vector<unsigned long > AdyacenciasVertex(unsigned long i); //Retorna las adyacencias de 1 nodo

  std::vector < std::vector < unsigned long > > Mapaavector();

  std::vector < std::vector < unsigned long > > FloydWarshall();
  std::deque< TValue > getFloyd(std::vector < std::vector < unsigned long > > r, unsigned long a, unsigned long b);

  void PrintPlaneGraph( ) const;
  void PrintPreorderGraph( unsigned long i, std::vector<bool>& m);
  void PrintLevelsGraph( unsigned long i ) const;
  void PrintGraphAsPNG( const std::string& filename ) const;

  bool isSafe(unsigned long v, unsigned long path[], unsigned long pos);
  bool hamCycleUtil(unsigned long path[], unsigned long pos);
protected:
  TVertices m_Vertices;
  TMatrix   m_Costs;
};

#include "Graph.hxx"

#endif // __Graph__h__

// eof - Graph.h

