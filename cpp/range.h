// A few helper methods to work with ranges (assumed to be STL containers: value_type, begin and end)
#pragma once

#include <vector>
#include <algorithm>

#include <boost/iterator/iterator_facade.hpp>

#include <ROOT/RVec.hxx>

namespace rdfhelpers {
template<typename RANGE,typename PREDICATE>
ROOT::VecOps::RVec<typename RANGE::value_type> select(const RANGE& range, PREDICATE&& pred)
{
  ROOT::VecOps::RVec<typename RANGE::value_type> sel;
  std::copy_if(range.begin(), range.end(), std::back_inserter(sel), pred);
  return sel;
}

template<typename RESULT,typename RANGE,typename FUNCTION>
typename ROOT::VecOps::RVec<RESULT> transform(const RANGE& range, FUNCTION&& fun)
{
  typename ROOT::VecOps::RVec<RESULT> result;
  std::transform(range.begin(), range.end(), std::back_inserter(result), fun);
  return result;
}

template<typename RANGE,typename PREDICATE>
typename RANGE::value_type next(const RANGE& range, PREDICATE&& pred)
{
  return *std::find_if(range.begin(), range.end(), pred);
}

template<typename RANGE,typename RESULT,typename REDUCE>
RESULT reduce(const RANGE& range, RESULT init, REDUCE&& reduce)
{
  return std::accumulate(range.begin(), range.end(), init, reduce);
}

/*
 * Provide an iterator interface over an index that runs from 0 to N
 */
template<typename IDX>
class IndexRangeIterator
  : public boost::iterator_facade<IndexRangeIterator<IDX>,
      IDX, // value
      boost::random_access_traversal_tag,
      IDX  // ref
    >
{
public:
  IndexRangeIterator(IDX idx, IDX max)
    : m_idx{idx}, m_max{max} {}

  IDX dereference() const { return m_idx; }
  bool equal( const IndexRangeIterator& other ) const { return ( m_max == other.m_max ) && ( m_idx == other.m_idx ); }
  void increment() { ++m_idx; }
  void decrement() { --m_idx; }
  void advance(std::ptrdiff_t d) { m_idx += d; }
  std::ptrdiff_t distance_to(const IndexRangeIterator& other) const { return m_idx-other.m_idx; }
private:
  IDX m_idx, m_max;
};
// fake container
template<typename IDX>
class IndexRange {
public:
  using iterator = IndexRangeIterator<IDX>;
  using const_iterator = IndexRangeIterator<IDX>;
  using value_type = IDX;

  explicit IndexRange(IDX max) : m_max(max) {}

  iterator begin() const { return IndexRangeIterator<IDX>(0, m_max); }
  iterator end  () const { return IndexRangeIterator<IDX>(m_max, m_max); }
private:
  IDX m_max;
};
};
