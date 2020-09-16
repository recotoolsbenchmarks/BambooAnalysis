// A few helper methods to work with ranges (assumed to be STL containers: value_type, begin and end)
#pragma once

#include <vector>
#include <algorithm>

#include <boost/iterator/iterator_facade.hpp>

#include <ROOT/RVec.hxx>
#include <TMath.h>

namespace rdfhelpers {

template<typename RANGE>
Double_t RMS(RANGE&& range)
{
    return TMath::RMS(range.begin(), range.end());
}


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
  if ( std::none_of(std::begin(lengths), std::end(lengths), [] ( std::size_t len ) { return len == 0; }) ) {
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
  IDX operator[] (std::size_t i) const { return m_min+i; }
private:
  IDX m_min;
  IDX m_max;
};

/*
 * Ranges and iterators for GenParticle and HepMC trees
 */
namespace gen {

template<typename PARENTIDXS>
class ancestor_iterator
  : public boost::iterator_facade<ancestor_iterator<PARENTIDXS>,
      int, // value
      boost::forward_traversal_tag,
      int // ref
    >
{
public:
  ancestor_iterator(int idx, const PARENTIDXS& parentIdxs)
    : m_idx(idx), m_parentIdxs(parentIdxs) {}

  int dereference() const { return m_idx; }
  bool equal( const ancestor_iterator& other ) const { return ( &m_parentIdxs == &other.m_parentIdxs ) && ( m_idx == other.m_idx ); }
  void increment() { m_idx = m_parentIdxs[m_idx]; }
private:
  int m_idx;
  const PARENTIDXS& m_parentIdxs;
};
template<int invalid,typename PARENTIDXS>
class ancestors
{
public:
  using iterator = ancestor_iterator<PARENTIDXS>;
  using const_iterator = ancestor_iterator<PARENTIDXS>;
  using value_type = int;

  ancestors(int self, const PARENTIDXS& parentIdxs)
    : m_self(self), m_parentIdxs(parentIdxs) {}

  iterator begin() const { return iterator(m_self , m_parentIdxs); }
  iterator end  () const { return iterator(invalid, m_parentIdxs); }
private:
  int m_self;
  const PARENTIDXS& m_parentIdxs;
};

template<int invalid,int end_offset,typename CHILDIDXS>
class descendant_iterator
  : public boost::iterator_facade<descendant_iterator<invalid,end_offset,CHILDIDXS>,
      int, // value
      boost::forward_traversal_tag,
      int // ref
    >
{
public:
  descendant_iterator(int self, const CHILDIDXS& childIdxs_begin, const CHILDIDXS& childIdxs_end)
    : m_idx(invalid), m_end(invalid), m_childIdxs_begin(childIdxs_begin), m_childIdxs_end(childIdxs_end) {
    if ( self != invalid ) {
      m_idx = m_childIdxs_begin[self];
      m_end = getEndChild(self);
    }
  }

  int dereference() const { return m_idx; }
  bool equal( const descendant_iterator& other ) const { return ( &m_childIdxs_begin == &other.m_childIdxs_begin ) && ( &m_childIdxs_end == &other.m_childIdxs_end ) && ( m_idx == other.m_idx ); }
  void increment() {
    // find descendants
    const auto c_begin = m_childIdxs_begin[m_idx];
    const auto c_end   = getEndChild(m_idx);
    if ( c_begin != c_end ) {
      if ( (m_idx+1) != m_end )
        m_stack_remain.emplace_back(m_idx+1, m_end);
      m_idx = c_begin;
      m_end = c_end;
    } else { // move to next sibling
      ++m_idx;
      if ( m_idx == m_end ) { // last -> move up
        if ( ! m_stack_remain.empty() ) {
          const auto nxt = m_stack_remain.back();
          m_idx = nxt.first;
          m_end = nxt.second;
          m_stack_remain.pop_back();
        } else { // end
          m_idx = invalid;
          m_end = invalid;
        }
      }
    }
  }
private:
  int m_idx, m_end; // current and end index (one past the last sibling)
  std::vector<std::pair<int,int>> m_stack_remain; // stack with remaining begin-end sibling ranges
  const CHILDIDXS& m_childIdxs_begin;
  const CHILDIDXS& m_childIdxs_end;

  int getEndChild(int parent) const {
    const auto end = m_childIdxs_end[parent];
    return end + ( ( end_offset && ( end != invalid ) ) ? end_offset : 0 );
  }
};
template<int invalid,int end_offset,typename CHILDIDXS>
class descendants_impl
{
public:
  using iterator = descendant_iterator<invalid,end_offset,CHILDIDXS>;
  using const_iterator = descendant_iterator<invalid,end_offset,CHILDIDXS>;
  using value_type = int;

  descendants_impl(int self, const CHILDIDXS& childIdxs_begin, const CHILDIDXS& childIdxs_end)
    : m_self(self), m_childIdxs_begin(childIdxs_begin), m_childIdxs_end(childIdxs_end) {}

  iterator begin() const { return iterator(m_self , m_childIdxs_begin, m_childIdxs_end); }
  iterator end  () const { return iterator(invalid, m_childIdxs_begin, m_childIdxs_end); }
private:
  int m_self;
  const CHILDIDXS& m_childIdxs_begin;
  const CHILDIDXS& m_childIdxs_end  ;
};
template<int invalid,typename CHILDIDXS> using descendants_firstlast = descendants_impl<invalid,1,CHILDIDXS>;
template<int invalid,typename CHILDIDXS> using descendants = descendants_impl<invalid,0,CHILDIDXS>;
};
};
