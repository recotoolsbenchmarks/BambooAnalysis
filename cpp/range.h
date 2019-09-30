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
  //std::copy_if(range.begin(), range.end(), std::back_inserter(sel), pred);
  // explicit copy to work around universal-reference push_back of RVec
  for ( const auto& itm : range ) {
    if ( pred(itm) ) {
      sel.push_back(typename RANGE::value_type(itm));
    }
  }
  return sel;
}

template<typename RANGE,typename PREDICATE>
ROOT::VecOps::RVec<typename RANGE::value_type> sort(const RANGE& range, PREDICATE&& pred)
{
  using item_t = typename RANGE::value_type;
  ROOT::VecOps::RVec<item_t> sel{std::begin(range), std::end(range)};
  std::stable_sort(std::begin(sel), std::end(sel),
      [&pred] ( const item_t ia, const item_t ib ) {
        return pred(ia) < pred(ib);
      });
  return sel;
}

template<typename RESULT,typename RANGE,typename FUNCTION>
typename ROOT::VecOps::RVec<RESULT> map(const RANGE& range, FUNCTION&& fun)
{
  typename ROOT::VecOps::RVec<RESULT> result;
  result.reserve(range.size());
  std::transform(range.begin(), range.end(), std::back_inserter(result), fun);
  return result;
}

template<typename RANGE,typename PREDICATE>
typename RANGE::value_type next(const RANGE& range, PREDICATE&& pred, const typename RANGE::value_type invalid)
{
  const auto it = std::find_if(range.begin(), range.end(), pred);
  return it != range.end() ? *it : invalid;
}

template<typename RANGE,typename RESULT,typename REDUCE>
RESULT reduce(const RANGE& range, RESULT init, REDUCE&& reduce)
{
  return std::accumulate(range.begin(), range.end(), init, reduce);
}

/**
 * N-particle combination
 *
 * TODO extend (could also carry common vertex and momentum)
 */
template<std::size_t NUM>
class Combination {
public:
  Combination() = default;
  Combination(const typename std::array<std::size_t,NUM> idx)
    : m_idx(idx) {}

  std::size_t get(std::size_t i) const { return m_idx.at(i); }
private:
  typename std::array<std::size_t,NUM> m_idx;
};

/**
 * 2-particle combinatorics from different base containers
 *
 * no duplicate or overlap checks are performed, these should be included in the predicate
 */
template<typename PREDICATE,typename RANGE1,typename RANGE2>
ROOT::VecOps::RVec<Combination<2>> combine2(PREDICATE&& pred, const RANGE1& range1, const RANGE2& range2)
{
  ROOT::VecOps::RVec<Combination<2>> sel;
  for ( const typename RANGE1::value_type i1 : range1 ) {
    for ( const typename RANGE2::value_type i2 : range2 ) {
      if ( pred(i1, i2) ) {
        sel.push_back(Combination<2>{std::array<std::size_t,2>{i1, i2}});
      }
    }
  }
  return sel;
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

  explicit IndexRange(IDX max) : m_min(  0), m_max(max) {}
  IndexRange(IDX min, IDX max) : m_min(min), m_max(max) {}

  iterator begin() const { return IndexRangeIterator<IDX>(m_min, m_max); }
  iterator end  () const { return IndexRangeIterator<IDX>(m_max, m_max); }
private:
  IDX m_min;
  IDX m_max;
};
};
