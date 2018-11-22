-- | The quaternions are a number system which generalise the complex numbers,
-- | in a similar way to how the complex numbers generalise the real numbers.
-- | One interesting aspect of the quaternions is that they satisfy all but one
-- | of the requirements to be a field; they only fail to be a field because
-- | their multiplication is non-commutative.
-- |
-- | Just as complex numbers can be thought of as pairs of real numbers, where
-- | the first is the "real part" and the second is the "imaginary part",
-- | quaternions can be thought of as a pair containing a "real part", which is
-- | again simply a real number, and a "vector part", which is a vector in R^3.
-- | Addition of quaternions is then easy to define: the sum of two quaternions
-- | has a real part equal to the sum of the operands' real parts, and a vector
-- | part equal to the sum of the operands' vector parts.
-- |
-- | There are four particularly important quaternions, which are known as the
-- | "basis elements". The first is `1`, which (when considered as a
-- | quaternion) has a real part of `1` and a vector part of `0` (i.e. the zero
-- | vector). The standard basis vectors in R^3 are commonly written as `e_1`,
-- | `e_2`, and `e_3`; the quaternions with a real part of 0 and a vector part
-- | of `e_1`, `e_2`, and `e_3` are called `i`, `j`, and `k` respectively, and
-- | form the other three basis elements.
-- |
-- | These four elements are called the basis elements because every quaternion
-- | can be written as `a + bi + cj + dk` for a choice of real numbers `a`,
-- | `b`, `c`, and `d`.
-- |
-- | Now, we can define multiplication on quaternions just by defining
-- | multiplication on the basis elements `1`, `i`, `j`, and `k`, for which we
-- | use the formula:
-- |
-- | ```
-- | i^2 = j^2 = k^2 = ijk = -1
-- | ```
-- |
-- | One of the first things we can deduce is that the multiplicative inverse
-- | of `i` is `-i` (just as with complex numbers), and also that the inverses
-- | of `j` and `k` are `-j` and `-k` respectively.
-- |
-- | From here we can calculate products between any two basis elements. For
-- | example, to calculate the product `jk`, we notice that `i^2 = ijk` and
-- | therefore we can cancel `i` on the left to yield `i = jk`. As another
-- | example, to calculate the product `ik`, we may start with the equation
-- | `i = jk` and multiply both sides by `k`, yielding `ik = jk^2 = -j`. All
-- | other products of basis elements may be obtained in a similar fashion.
-- |
-- | Knowing the products of basis elements and using the distributive law, we
-- | can find the product of any two quaternions:
-- |
-- | ```
-- | (a + bi + cj + dk) * (e + fi + gj + hk)
-- | = ae + afi + agj + ahk
-- |   + bei + bf(i^2) + bg(ij) + bh(ik)
-- |   + cej + cf(ji) + cg(j^2) + ch(jk)
-- |   + dek + df(ki) + dg(kj) + dh(k^2)
-- | = ae - bf - cg - dh
-- |   + (af + be + ch - dg) i
-- |   + (ag - bh + ce + df) j
-- |   + (ah + bg - cf + de) k
-- | ```
-- |
-- | Note that quaternion multiplication is not commutative; that is, `p * q`
-- | is usually not the same as `q * p`. For example, `ij = k`, but `ji = -k`.
-- |
-- | Like the real numbers, however, each quaternion does have a multiplicative
-- | inverse, i.e. for each quaternion `p` there exists a unique quaternion `q`
-- | such that `p * q = q * p = 1`.
-- |
-- | This means that there are two ways of dividing quaternions. If we want to
-- | divide a quaternion `p` by another quaternion `q`, we can multiply `p` by
-- | the multiplicative inverse of `q`, which can be written `q^-1`. Of course,
-- | we have two ways of doing this; `p * q^-1` and `q^-1 * p` will not
-- | necessarily be the same, so we are really dealing with two different
-- | operations here. We call the operation taking `p` and `q` to `q^-1 * p`
-- | "left-division", and the alternative operation which yields `p * q^-1`
-- | "right-division".
-- |
-- | The most important application of quaternions in computing is for
-- | representing orientations and rotations in 3D space; see the
-- | Data.Quaternion.Rotation module for more details of this.
module Data.Quaternion
  ( module Data.Quaternion
  , module ReExports
  ) where

import Prelude
import Math as Math
import Data.Quaternion.Vec3 (Vec3, vec3)

import Data.DivisionRing (leftDiv, rightDiv) as ReExports

-- | A quaternion. The type parameter denotes the underlying type. Note that
-- | the underlying type should be a reasonable approximation of the real
-- | numbers; if this is not the case, some of the functions may exhibit
-- | strange behaviour.
-- |
-- | Because multiplication of quaternions is non-commutative, there is no
-- | `CommutativeRing` instance, and consequently no `EuclideanRing` or `Field`
-- | instance either.
-- |
-- | This means, amongst other things, that the `(/)` operator from Prelude
-- | cannot be used with `Quaternion` values. However, `Quaternion` does have
-- | a `DivisionRing` instance, so you can use `leftDiv` and `rightDiv` from
-- | the module `Data.DivisionRing` from the `prelude` library instead. These
-- | functions are also re-exported from this module for convenience.
data Quaternion a = Quaternion a a a a

derive instance eqQuaternion :: Eq a => Eq (Quaternion a)

instance showQuaternion :: Show a => Show (Quaternion a) where
  show (Quaternion a b c d) =
    "(Quaternion " <>
      show a <> " " <>
      show b <> " " <>
      show c <> " " <>
      show d <> ")"

instance semiringQuaternion :: Ring a => Semiring (Quaternion a) where
  zero =
    Quaternion zero zero zero zero
  add (Quaternion a1 b1 c1 d1) (Quaternion a2 b2 c2 d2) =
    Quaternion (a1 + a2) (b1 + b2) (c1 + c2) (d1 + d2)
  one =
    Quaternion one zero zero zero
  mul (Quaternion a1 b1 c1 d1) (Quaternion a2 b2 c2 d2) =
    Quaternion
      (a1*a2 - b1*b2 - c1*c2 - d1*d2)
      (a1*b2 + b1*a2 + c1*d2 - d1*c2)
      (a1*c2 - b1*d2 + c1*a2 + d1*b2)
      (a1*d2 + b1*c2 - c1*b2 + d1*a2)

instance ringQuaternion :: Ring a => Ring (Quaternion a) where
  sub (Quaternion a1 b1 c1 d1) (Quaternion a2 b2 c2 d2) =
    Quaternion (a1 - a2) (b1 - b2) (c1 - c2) (d1 - d2)


instance functorQuaternion :: Functor Quaternion where
  map f (Quaternion a b c d) = Quaternion (f a) (f b) (f c) (f d)

instance divisionRingQuaternion :: DivisionRing a => DivisionRing (Quaternion a) where
  recip q = scalarMul (recip (normSquare q)) (conjugate q)


-- | The real part of the quaternion, that is, the first component. Defined as
-- |
-- |     \(Quaternion a _ _ _) -> a
-- |
realPart :: forall a. Quaternion a -> a
realPart (Quaternion a _ _ _) = a

-- | The vector part of the quaternion, that is, the second, third, and fourth
-- | components, represented as an array with exactly 3 elements. Defined as
-- |
-- |     \(Quaternion _ x y z) -> vec3 x y z
-- |
vectorPart :: forall a. Quaternion a -> Vec3 a
vectorPart (Quaternion _ x y z) = vec3 x y z

-- | The conjugate of a quaternion. This operation negates the vector part of
-- | the quaternion.
conjugate :: forall a. Ring a => Quaternion a -> Quaternion a
conjugate (Quaternion a b c d) =
  Quaternion a (-b) (-c) (-d)

-- | The conjugate of a quaternion by another quaternion. Defined as
-- |
-- |     \p q -> q * p * recip q
-- |
conjugateBy :: forall a. DivisionRing a => Quaternion a -> Quaternion a -> Quaternion a
conjugateBy p q = q * p * recip q

norm :: Quaternion Number -> Number
norm q = Math.sqrt (normSquare q)

-- | The square of the norm of a quaternion. This is slightly easier to compute
-- | than the actual norm, so may be useful in cases where you are worried
-- | about performance.
normSquare :: forall a. Semiring a => Quaternion a -> a
normSquare (Quaternion a b c d) = a*a + b*b + c*c + d*d

i :: forall a. Semiring a => Quaternion a
i = Quaternion zero one zero zero

j :: forall a. Semiring a => Quaternion a
j = Quaternion zero zero one zero

k :: forall a. Semiring a => Quaternion a
k = Quaternion zero zero zero one

-- | Multiplies both the real part and the vector part by the given scalar.
scalarMul :: forall a. Semiring a => a -> Quaternion a -> Quaternion a
scalarMul k' (Quaternion a b c d) =
  Quaternion (k' * a) (k' * b) (k' * c) (k' * d)

-- | Scales the given quaternion, returning a quaternion pointing in the same
-- | direction of unit norm; multiplying this by the original quaternion's
-- | norm will give you back the original quaternion.
versor :: Quaternion Number -> Quaternion Number
versor q = scalarMul (1.0 / norm q) q

-- | Approximate equality of quaternions, given an epsilon value specifying the
-- | maximum amount that any of the four components is allowed to differ by.
approxEq :: forall a. Ord a => Ring a => a -> Quaternion a -> Quaternion a -> Boolean
approxEq eps (Quaternion a1 b1 c1 d1) (Quaternion a2 b2 c2 d2) =
  ok a1 a2 && ok b1 b2 && ok c1 c2 && ok d1 d2
  where
  ok x y =
    let
      z = x - y
      diff = if z > zero then z else negate z
    in
      diff < eps
