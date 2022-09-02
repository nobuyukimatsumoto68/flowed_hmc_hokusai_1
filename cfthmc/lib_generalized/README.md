## Flow kernels made with Wilson loops

Applicable flow type (9/21 present):

[1] Kernel: made from Wilson loops, whose dependence on the link variables is linear. We divide a flow step into several pieces and make updates with them separately. In particular, links directed in different directions are updated separately. We say that each piece of the flow step has different "flow type". It is assumed below that Wilson loops has the same length for each flow type.

[2] Masking scheme: arbitrary (assuming that the whole lattice can be filled with the same pattern). Different masking can be used for different flow types.

### Strategy

- We construct a table of involved links for each flow type. We make use of the following symmetries:
  - rotation: storing the table for $\mu$ (flowed direction) $=0$ case suffices.
  - translation: flowed link can be set in the fundamental unit including (0,0,0,0).
- We make use of the table when calculating derivatives of Wilson loops. In particular, the second derivatives will be required in evaluating Jacobi matrices. 

### Terminologies used in the code

 **(fundamental) unit**: The unit defined in making a specific symmetry pattern in [2]; this will be copied (without rotations) to fill the whole lattice. 

​	e.g., ordinary even/odd masking:
$$
{e/o} = (x_1+x_2+x_3+x_4)\%2
$$
​	fundamental unit: 2x2x2x2 lattice

**index**: an integer which labels the sites in the unit.

**block**: The coordinates of the units. These blocks are parameterized by four integers.

**shape**: vector of integers representing the shape of a Wilson loop. It is assumed that the Wilson loops have the same shapes for all sites. We do not store inversed loops assuming hermiticity.

### user-defined classes

- **UnitInfo**: stores information of the fundamental unit. 
- **KerneInfo**: stores information of kernels. 
- **MaskInfo**: stores information of masking scheme. 
- **FlowType**: stores the table. 
  - each row of the table has a structure defined in **RowOfTable**

### Example: plaquette flow with e/o masking

- **Declare the flow type:** 

  ```c++
  int mu=0; // flowed direction
  std::array<int,4> unit_dimensions = std::array<int,4>{{2,2,2,2}}; // unit
  int n_loops = 6; // number of loops (shapes)
  int length_of_loop = 4; // length of loops
  int n_masks = 2; // number of masking types
  int description = 1; // distinct integer which we assign to each flow type
  int flow_size = 1; // size to be extended for MPI communications
  
  FlowType plaq(mu, unit_dimensions, n_loops, length_of_loops, n_masks,
                description, flow_size);
  ```

- **Set shapes**

  ```c++
  for(int nu=mu+1; nu<4; ++nu)
      plaq.set_shape(nu-1,std::vector<int>{mu,nu,-mu-1,-nu-1});
  for(int nu=mu+1; nu<4; ++nu)
      plaq.set_shape(3+nu-1,std::vector<int>{mu,-nu-1,-mu-1,nu});
  ```

- **Set mask pattern**


  ```c++
  const std::function<int(const int, const int, const int, const int)> even_odd
      = [](const int x1, const int x2, const int x3, const int x4) {
      return (x1 + x2 + x3 + x4 )%2; };
  plaq.set_mask( even_odd );
  ```

- **Calculate table + supply the flow type to std::unordered_map (for a technical reason)**

  ```c++
  plaq.initialize();
  GlobalScopeHash::get_flow_type.insert({ plaq.m_description, plaq });
  ```

This is done in 

```c++
GlobalScopeHash::setup_hash_map();
```



In the plaquette case, we obtain the following table:

```
x: idx =  0, mask =  0, y: idx =  0, nu =  0, mask =  0, block =  0,  0,  0,  0, is_positive =  1, l =  5, s =  0, m_prime =  0, next =  6
x: idx =  0, mask =  0, y: idx =  0, nu =  0, mask =  0, block =  0,  0,  0,  0, is_positive =  1, l =  1, s =  0, m_prime =  0, next =  5
x: idx =  0, mask =  0, y: idx =  0, nu =  0, mask =  0, block =  0,  0,  0,  0, is_positive =  1, l =  0, s =  0, m_prime =  0, next =  4
x: idx =  0, mask =  0, y: idx =  0, nu =  0, mask =  0, block =  0,  0,  0,  0, is_positive =  1, l =  2, s =  0, m_prime =  0, next =  3
x: idx =  0, mask =  0, y: idx =  0, nu =  0, mask =  0, block =  0,  0,  0,  0, is_positive =  1, l =  4, s =  0, m_prime =  0, next =  2
x: idx =  0, mask =  0, y: idx =  0, nu =  0, mask =  0, block =  0,  0,  0,  0, is_positive =  1, l =  3, s =  0, m_prime =  0, next =  1
x: idx =  0, mask =  0, y: idx =  0, nu =  1, mask =  0, block =  0,  0,  0,  0, is_positive =  0, l =  0, s =  3, m_prime =  1, next =  1
x: idx = 12, mask =  0, y: idx =  0, nu =  1, mask =  0, block =  1,  0,  0,  0, is_positive =  0, l =  3, s =  1, m_prime =  2, next =  1
x: idx =  0, mask =  0, y: idx =  0, nu =  2, mask =  0, block =  0,  0,  0,  0, is_positive =  0, l =  1, s =  3, m_prime =  3, next =  1
x: idx = 10, mask =  0, y: idx =  0, nu =  2, mask =  0, block =  1,  0,  0,  0, is_positive =  0, l =  4, s =  1, m_prime =  4, next =  1
x: idx =  0, mask =  0, y: idx =  0, nu =  3, mask =  0, block =  0,  0,  0,  0, is_positive =  0, l =  2, s =  3, m_prime =  5, next =  1
x: idx =  9, mask =  0, y: idx =  0, nu =  3, mask =  0, block =  1,  0,  0,  0, is_positive =  0, l =  5, s =  1, m_prime =  6, next =  1
x: idx =  0, mask =  0, y: idx =  1, nu =  0, mask =  1, block =  0,  0,  0, -1, is_positive =  0, l =  5, s =  2, m_prime =  0, next =  1
x: idx =  0, mask =  0, y: idx =  1, nu =  0, mask =  1, block =  0,  0,  0,  0, is_positive =  0, l =  2, s =  2, m_prime =  1, next =  1
x: idx =  3, mask =  0, y: idx =  1, nu =  0, mask =  1, block =  0,  0,  0,  0, is_positive =  0, l =  4, s =  2, m_prime =  2, next =  1
x: idx =  3, mask =  0, y: idx =  1, nu =  0, mask =  1, block =  0,  0,  1,  0, is_positive =  0, l =  1, s =  2, m_prime =  3, next =  1
x: idx =  5, mask =  0, y: idx =  1, nu =  0, mask =  1, block =  0,  0,  0,  0, is_positive =  0, l =  3, s =  2, m_prime =  4, next =  1
x: idx =  5, mask =  0, y: idx =  1, nu =  0, mask =  1, block =  0,  1,  0,  0, is_positive =  0, l =  0, s =  2, m_prime =  5, next =  1
...
```

$(x, \mu: fixed)$ represents the flowed link and $(y, \nu)$ the link which is involved in the flow. "is_pos" shows if the link variable on $(y, \nu)$ is positive directed when that on $(x, \mu)$ is positive-directed. $\ell$ labels the shapes and $s$ labels the sites in a given loop $\ell$. $m'$ labels distinct yidx for a given xidx.

The map between coordinates and x_idx are given by

```
x_unit = (0, 0, 0, 0). x_idx = 0
x_unit = (0, 0, 0, 1). x_idx = 1
x_unit = (0, 0, 1, 0). x_idx = 2
x_unit = (0, 0, 1, 1). x_idx = 3
x_unit = (0, 1, 0, 0). x_idx = 4
x_unit = (0, 1, 0, 1). x_idx = 5
x_unit = (0, 1, 1, 0). x_idx = 6
x_unit = (0, 1, 1, 1). x_idx = 7
x_unit = (1, 0, 0, 0). x_idx = 8
x_unit = (1, 0, 0, 1). x_idx = 9
x_unit = (1, 0, 1, 0). x_idx = 10
x_unit = (1, 0, 1, 1). x_idx = 11
x_unit = (1, 1, 0, 0). x_idx = 12
x_unit = (1, 1, 0, 1). x_idx = 13
x_unit = (1, 1, 1, 0). x_idx = 14
x_unit = (1, 1, 1, 1). x_idx = 15
```

The table is sorted by the following rule:

```c++
class RowOfTable{
    static constexpr int kXIdx= 0;
	static constexpr int kXMask= 1;
	static constexpr int kYIdx= 2;
	static constexpr int kNu= 3;
	static constexpr int kYMask= 4;
	static constexpr int kYBlockBeg= 5;
  	// 5--8 : block
	static constexpr int kIsPositive= 9;
	static constexpr int kL= 10;
	static constexpr int kS= 11;
	static constexpr int kMPrime= 12;
    static constexpr int kNext= 13;
...
}

bool RowOfTable::operator< (const RowOfTable& another) const {

  int i=kXMask; // x_mask
  if( (*this)[i] < another[i] ) return true;
  else if( (*this)[i] > another[i] ) return false;

  i=kYIdx; // y_idx
  if( (*this)[i] < another[i] ) return true;
  else if( (*this)[i] > another[i] ) return false;

  i=kNu; // nu
  if( (*this)[i] < another[i] ) return true;
  else if( (*this)[i] > another[i] ) return false;

  i=kXIdx; // x_idx
  if( (*this)[i] < another[i] ) return true;
  else if( (*this)[i] > another[i] ) return false;

  for(int rho=0; rho<4; ++rho){
    i=kYBlockBeg+rho; // y_blocks
    if( (*this)[i] < another[i] ) return true;
    else if( (*this)[i] > another[i] ) return false;
  }

  i=kIsPositive; // is_positive
  if( (*this)[i] < another[i] ) return true;
  else if( (*this)[i] > another[i] ) return false;

  i=kYMask; // y_mask
  if( (*this)[i] < another[i] ) return true;
  else if( (*this)[i] > another[i] ) return false;

  i=kMPrime; // m_prime
  if( (*this)[i] < another[i] ) return true;
  else if( (*this)[i] > another[i] ) return false;
    
  i=kNext; // next
  if( (*this)[i] < another[i] ) return true;
  else if( (*this)[i] > another[i] ) return false;

  return false;
}
```

We iterate over the pairs $(x, \nu)$ for a given $(y, \mu: fixed)$ in the code. With the ordering above, this can be done by iterating over the rows of the table.

The first row with the given (xmask, yidx) can be obtained by:

```c++
inline int begin_mask_yidx(const int mask, const int y_idx);
```



For a general flow direction $\mu$, we first need to make a rotation. This can be done by calling

```c++
// coordinates
inline std::array<int, 4> rotate( const std::array<int, 4>& x, const int mu ) const;
// direction
inline int rotate( const int nu, const int mu) const;
```

e.g.,

```c++
const Coordinate yg = ...;
const std::array<int, 4> y_prime = flow_type.rotate( Coordinate2std(yg), -mu-1 ); 
```

After referring to the table, we return the direction back:

```c++
const int x_prime_idx = flow_type.table(pos,0);
const int nu_prime = flow_type.table(pos,3);

const Coordinate xg = std2Coordinate( flow_type.rotate( x_prime, mu ) );
const int nu = flow_type.rotate( nu_prime, mu ); 
```

