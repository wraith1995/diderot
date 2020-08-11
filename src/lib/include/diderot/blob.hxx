/*! \file tensor.hxx
 *
 * \author Teo Collin
 *
 * Base classes for the generated cache types in the compiler
 * 
 */

/*
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 */

#ifndef _DIDEROT_BLOB_HXX_
#define _DIDEROT_BLOB_HXX_

namespace diderot {

  template <typename REAL, const int N>
  struct blob_real {
    REAL _data[N];
    blob_real ()
    {
      REAL a = 0.0;
      std::fill_n(this->_data, N, a);
    }
    REAL * copyout (REAL * data, const int start, const int size){
      std::memcpy(data, & (_data[start]), size * sizeof (REAL));
      return data;
    }
    void copyin (const REAL * data, const int start, const int size){
      std::memcpy(& (_data[start]), data, size * sizeof (REAL));
      return;
    }
    inline REAL getscalar (const int start, const int offset){
      return _data[start + offset];
    }
  };
  template <typename REAL, const int N>
  struct blob_int {
    REAL _data[N];
    blob_int ()
    {
      REAL a = -1;
      std::fill_n(this->_data, N, a);
    }
    bool check(const int addr, REAL val)
    {
      return (_data[addr] == val);
    }
    void set(const int addr, REAL val){
      _data[addr] = val;
      return;
    }
  };

} // namespace diderot

#endif // !_DIDEROT_TENSOR_HXX_
