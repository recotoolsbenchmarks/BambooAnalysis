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

namespace detail {
template<std::size_t N,typename PRED,std::size_t... IDX,typename... RANGES>
void combine_add_if(ROOT::VecOps::RVec<Combination<N>>& out, PRED&& pred, const std::array<std::size_t,N>& indices, std::index_sequence<IDX...>, RANGES&&... ranges)
{
  if ( pred(ranges[indices[IDX]]...) ) {
    out.push_back(Combination<N>({ ranges[indices[IDX]]... }));
  }
}
}
/**
 * N-particle combinatorics
 *
 * no duplicate or overlap checks are performed, these should be included in the predicate
 */
template<typename PREDICATE,typename... RANGES>
ROOT::VecOps::RVec<Combination<sizeof...(RANGES)>> combine(PREDICATE&& pred, RANGES&&... ranges)
{
  constexpr auto N = sizeof...(RANGES);
  constexpr auto idx_seq = std::make_index_sequence<N>{}; // to zip idx and ranges
  ROOT::VecOps::RVec<Combination<N>> out;
  const auto lengths = std::array<size_t,N>{{ranges.size()...}};
  std::array<std::size_t,N> idx{};
  idx.fill(0);
  std::size_t j = N; // index+1 of the outermost array that is last updated
  while ( j != 0 ) {
    detail::combine_add_if(out, pred, idx, idx_seq, ranges...);
    // increase the (N-dimensional) counter (starting from the last array, as a nested loop would do)
    j = N;
    while ( j > 0 ) {
      if ( idx[j-1]+1 != lengths[j-1] ) {
        ++idx[j-1];
        break;
      } else
        idx[--j] = 0;
    }
  }
  return out;
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
  std::size_t size() const { return m_max-m_min; }
private:
  IDX m_min;
  IDX m_max;
};
};
