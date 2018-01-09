/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012   FSU Statistical Shape Analysis and Modeling Group
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#ifndef KDTREE_H
#define KDTREE_H 1

#include <cstddef>
#include <limits>

namespace srvf
{

/**
 * k-d tree implementation.
 *
 * This should work for any datatype \c T which supports \c operator[], 
 * \c size(), and \c distance_to().  Currently, it is used to store 
 * \c Matrix objects for the partial matching routines.
 */
template <class T>
class KdTree
{
private:
  
  struct Node
  {
    Node ()
     : parent(NULL), left(NULL), right(NULL)
    { }

    Node (const T &val, Node *parent, Node *left, Node *right)
     : data(val), parent(parent), left(left), right(right)
    { }

    T     data;
    Node *parent;
    Node *left;
    Node *right;
  };

public:

  /**
   * Creates an empty k-d tree.
   */
  KdTree()
   : root_(NULL)
  { }

  ~KdTree()
  {
    kdt_free_(root_);
  }

  /**
   * Inserts an element into the tree.
   */
  void insert(const T &newval)
  {
    if (root_ == NULL)
      root_ = new Node(newval, NULL, NULL, NULL);
    else
      kdt_insert_(newval, root_, 0);
  }

  /**
   * Inserts the given element into the tree if and only if it is at least 
   * \a min_dist away from all elements in the tree.
   *
   * \return \c true if a new element was inserted, or \c false otherwise
   */
  bool insert_cond(const T &newval, double min_dist)
  {
    if (root_ == NULL){
      root_ = new Node(newval, NULL, NULL, NULL);
      return true;
    }
    else
    {
      T best_nbr;
      double best_dist=1e9;
      kdt_get_nn_(newval, root_, 0, &best_nbr, &best_dist);
      if (best_dist >= min_dist)
      {
        kdt_insert_(newval, root_, 0);
        return true;
      }
      else
      {
        return false;
      }
    }
  }

  /**
   * Finds the nearest neighbor of the given element.
   * 
   * If \a dist_ptr is non-null, it will receive the distance of 
   * the nearest neighbor.
   */
  T get_nearest_neighbor(const T &val, double *dist_ptr=NULL)
  {
    T best_nbr;
    double best_dist=std::numeric_limits<double>::infinity();

    kdt_get_nn_ (val, root_, 0, &best_nbr, &best_dist );
    if (dist_ptr != NULL) *dist_ptr = best_dist;
    return best_nbr;
  }

  std::vector<T> to_vector()
  {
    std::vector<T> result;
    kdt_to_vector_(root_, result);
    return result;
  }


private:

  void kdt_insert_(const T &newval, Node *root, size_t level)
  {
    size_t dim = level % newval.size();
    double x1 = newval[dim];
    double x2 = (root->data)[dim];

    if ( x1 <= x2 ){
      if ( root->left != NULL ){
        kdt_insert_( newval, root->left, level+1 );
      } else {
        root->left = new Node(newval, root, NULL, NULL);
      }
    } else {
      if ( root->right != NULL ){
        kdt_insert_( newval, root->right, level+1 );
      } else {
        root->right = new Node(newval, root, NULL, NULL);
      }
    }
  }

  void kdt_get_nn_ ( const T &A, Node *root, size_t level, 
    T *best_nbr, double *best_dist )
  {
    if ( !root ) return;

    double d = A.distance_to(root->data);
    if ( d < *best_dist )
    {
      *best_nbr = root->data;
      *best_dist = d;
    }

    size_t dim = level % A.size();
    double x1 = A[dim];
    double x2 = (root->data)[dim];

    Node *near_node, *far_node;
    if ( x1 <= x2 ){ near_node=root->left; far_node=root->right; }
    else           { near_node=root->right; far_node=root->left; }

    kdt_get_nn_( A, near_node, level+1, best_nbr, best_dist );

    if ( fabs( x1-x2 ) < *best_dist )
      kdt_get_nn_( A, far_node, level+1, best_nbr, best_dist );
  }

  void kdt_to_vector_(Node *root, std::vector<T> &v)
  {
    if (root != NULL)
    {
      kdt_to_vector_(root->left, v);
      v.push_back(root->data);
      kdt_to_vector_(root->right, v);
    }
  }

  void kdt_free_( Node *root ){
    if ( root != NULL )
    {
      kdt_free_( root->left );
      kdt_free_( root->right );
      delete root;
    }
  }

  Node *root_;
};

} // namespace srvf

#endif // KDTREE_H
