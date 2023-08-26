# PGA3D
3D Projective Geometric Algebra in C++.

## Install & use

To run the example:

```sh
g++ demo.cpp PGA3D.cpp -o demo
./demo
```

To run the tests:

```sh
g++ test.cpp PGA3D.cpp -o test
./test
```


## Implementation

The main actor is the `PGA3D` class, whose objects are called vectors.

Each vector consists of an array `mvec` of 16 floating-point numbers.
An `mvec` is a generalisation of the ordinary 3-vectors of linear algebra, as an `mvec` can also represent lines, planes, directions etc.

An `mvec` is solely a method of storing data, and in itself has no geometrical meaning.
The meaning is derived from the projective geometric algebra that is defined on them, i.e. a set of functions that defines operations between `mvec`s.

## Mathematical background

This section is mainly for my own quick reference.

### Algebra

- Geometric product: $a b$
- Dual: $a^*$
- Outer product (meet): $a \wedge b$
- Regressive product (join): $a \vee b$
- Inner product: $a \cdot b$
- Commutator product: $a \times b$

### Geometry

Geomtric primitives:

- Point $(x,y,z)$: $P = x e_{032} + y e_{013} + z e_{021} + e_{123}$
- Direction $(x,y,z)$: $D = x e_{032} + y e_{013} + z e_{021}$
- Plane $ax + by + cz + d = 0$: $P = a e_{1} + b e_{2} + c e_{3} + d e_{0}$

Making points:

- Meet planes $p_1$, $p_2$, $p_3$ in point $P$: $P = p_1 \wedge p_2 \wedge p_3$
- Meet line $\ell$ and plane $p$ in point $P$: $P = \ell \wedge p$

Making lines:

- Line $\ell$ through points $P_1$ and $P_2$: $\ell = P_1 \vee P_2$
- Meet planes $p_1$ and $p_2$ in line $\ell$: $\ell = p_1 \wedge p_2$

## Acknowledgements

- Code based on output of the [code generator for PGA3D](https://bivector.net/tools.html?p=3&q=0&r=1) on [bivector.net](https://bivector.net).
- Mathematical background from the [3DPGA cheatsheet](https://bivector.net/3DPGA.pdf) by C. Gunn, S. De Keninck

## License

This work is licensed under the MIT license. See [LICENSE](LICENSE) for details.
