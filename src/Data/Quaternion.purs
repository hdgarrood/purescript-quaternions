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
  ( Quaternion(..)
  , realPart
  , fromReal
  , vectorPart
  , fromVector
  , i
  , j
  , k
  , scalarMul
  , dot
  , conjugate
  , conjugateBy
  , approxEq
  , exp
  , log
  , norm
  , normSquare
  , infinityNorm
  , versor
  , module ReExports
  ) where

import Prelude

import Data.DivisionRing (leftDiv, rightDiv) as ReExports
import Data.Foldable (class Foldable)
import Data.Newtype (unwrap)
import Data.Ord (abs)
import Data.Ord.Max (Max(..))
import Data.Quaternion.Vec3 (Vec3, vec3)
import Data.Quaternion.Vec3 as Vec3
import Data.Semigroup.Foldable (class Foldable1, foldMap1)
import Math as Math
import Partial.Unsafe (unsafePartial)

-- | A quaternion. The type parameter denotes the underlying type. Note that
-- | the underlying type should be a reasonable approximation of the real
-- | numbers; if this is not the case, some of the functions may exhibit
-- | strange behaviour.
-- |
-- | The first argument to the `Quaternion` constructor is the real part, and
-- | the other three components constitute the vector part.
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
  show (Quaternion w x y z) =
    "(Quaternion " <>
      show w <> " " <>
      show x <> " " <>
      show y <> " " <>
      show z <> ")"

instance functorQuaternion :: Functor Quaternion where
  map f (Quaternion w x y z) = Quaternion (f w) (f x) (f y) (f z)

-- | The `Foldable` instance for `Quaternion` behaves as if the given
-- | `Quaternion` were converted to a four-element array before folding, like
-- | so:
-- |
-- |     foldr f z (Quaternion w x y z) = foldr f z [w, x, y, z]
-- |
-- | For example:
-- |
-- |     foldr (+) 0 (Quaternion 1.0 1.0 1.0 1.0) = 4.0
-- |
instance foldableQuaternion :: Foldable Quaternion where
  foldr f z' (Quaternion w x y z) =
    w `f` (x `f` (y `f` (z `f` z')))
  foldl f z' (Quaternion w x y z) =
    (((z' `f` w) `f` x) `f` y) `f` z
  foldMap = foldMap1

-- | The `Foldable1` instance for `Quaternion` behaves as if the given
-- | `Quaternion` were converted to a four-element array before folding. For
-- | example:
-- |
-- |     foldMap1 Additive (Quaternion 1.0 1.0 1.0 1.0) = Additive 4.0
-- |
instance foldable1Quaternion :: Foldable1 Quaternion where
  fold1 (Quaternion w x y z) = w <> x <> y <> z
  foldMap1 f (Quaternion w x y z) = f w <> f x <> f y <> f z

instance semiringQuaternion :: Ring a => Semiring (Quaternion a) where
  zero =
    Quaternion zero zero zero zero
  add (Quaternion w1 x1 y1 z1) (Quaternion w2 x2 y2 z2) =
    Quaternion (w1 + w2) (x1 + x2) (y1 + y2) (z1 + z2)
  one =
    Quaternion one zero zero zero
  mul (Quaternion w1 x1 y1 z1) (Quaternion w2 x2 y2 z2) =
    Quaternion
      (w1*w2 - x1*x2 - y1*y2 - z1*z2)
      (w1*x2 + x1*w2 + y1*z2 - z1*y2)
      (w1*y2 - x1*z2 + y1*w2 + z1*x2)
      (w1*z2 + x1*y2 - y1*x2 + z1*w2)

instance ringQuaternion :: Ring a => Ring (Quaternion a) where
  sub (Quaternion w1 x1 y1 z1) (Quaternion w2 x2 y2 z2) =
    Quaternion (w1 - w2) (x1 - x2) (y1 - y2) (z1 - z2)

instance divisionRingQuaternion :: DivisionRing a => DivisionRing (Quaternion a) where
  recip q = scalarMul (recip (normSquare q)) (conjugate q)

-- | The real part of the quaternion, that is, the first component. Defined as
-- |
-- |     \(Quaternion w _ _ _) -> w
-- |
realPart :: forall a. Quaternion a -> a
realPart (Quaternion w _ _ _) = w

-- | Construct a real quaternion (that is, a quaternion with vector part equal
-- | to the zero vector) from a single real number. This gives us the identity
-- |
-- |     fromReal (realPart q) + fromVector (vectorPart q) == q
-- |
fromReal :: forall a. Semiring a => a -> Quaternion a
fromReal w = Quaternion w zero zero zero

-- | The vector part of the quaternion, that is, the second, third, and fourth
-- | components, represented as an array with exactly 3 elements. Defined as
-- |
-- |     \(Quaternion _ x y z) -> vec3 x y z
-- |
vectorPart :: forall a. Quaternion a -> Vec3 a
vectorPart (Quaternion _ x y z) = vec3 x y z

-- | Construct an imaginary quaternion (that is, a quaternion with real part
-- | equal to zero) from a vector. The returned quaternion's vector part will
-- | be equal to the argument supplied. This gives us the identity
-- |
-- |     fromReal (realPart q) + fromVector (vectorPart q) == q
-- |
fromVector :: forall a. Semiring a => Vec3 a -> Quaternion a
fromVector v =
  unsafePartial
    case Vec3.toArray v of
      [x, y, z] -> Quaternion zero x y z

-- | The conjugate of a quaternion. This operation negates the vector part of
-- | the quaternion. Defined as
-- |
-- |     \(Quaternion w x y z) ->
-- |       Quaternion w (-x) (-y) (-z)
conjugate :: forall a. Ring a => Quaternion a -> Quaternion a
conjugate (Quaternion w x y z) =
  Quaternion w (-x) (-y) (-z)

-- | The conjugate of a quaternion by another quaternion. Defined as
-- |
-- |     \p q -> q * p * recip q
-- |
conjugateBy :: forall a. DivisionRing a => Quaternion a -> Quaternion a -> Quaternion a
conjugateBy p q = q * p * recip q

-- | The dot product of quaternions. Defined as
-- |
-- |     \(Quaternion w1 x1 y1 z1) (Quaternion w2 x2 y2 z2) ->
-- |       w1*w2 + x1*x2 + y1*y2 + z1*z2
dot :: forall a. Semiring a => Quaternion a -> Quaternion a -> a
dot (Quaternion w1 x1 y1 z1) (Quaternion w2 x2 y2 z2) =
  w1*w2 + x1*x2 + y1*y2 + z1*z2

-- | The quaternion exponential function. This function can be defined via a
-- | power series in a similar way to the real or complex exponential
-- | functions:
-- |
-- |     exp q = 1 + q + (q^2/2!) + (q^3/3!) + ...
-- |
-- | and if we let `q = a + v`, where `a` is the real part of `q` and `v` is
-- | the vector part of `q`, then for nonzero `v`, it can be shown that this
-- | definition can also be expressed as:
-- |
-- |     exp q = exp (a + v) = exp(a) * (cos (norm v) + (v/norm v) * sin (norm v))
-- |
-- | (see <https://math.stackexchange.com/questions/1030737/exponential-function-of-quaternion-derivation>).
-- | The implementation uses the second definition.
-- |
-- | For real quaternions &mdash; that is, quaternions with a vector part equal
-- | to the zero vector &mdash; this function coincides with the standard real
-- | exponential function. More generally, for complex quaternions &mdash; that
-- | is, quaternions with y and z component equal to zero &mdash; this function
-- | coincides with the complex exponential function.
-- |
-- | Note that due to noncommutativity of quaternion multiplication, the
-- | expected identity `exp (p + q) = exp p * exp q` will not necessarily hold
-- | unless `p` and `q` commute.
exp :: Quaternion Number -> Quaternion Number
exp q@(Quaternion w x y z) =
  let
    v = Quaternion 0.0 x y z
    normV = norm v
    k' = if normV == 0.0 then 0.0 else Math.sin normV / normV
  in
    scalarMul (Math.exp w)
      (Quaternion (Math.cos normV) (k'*x) (k'*y) (k'*z))

-- | The quaternion logarithm function. This function is the right inverse of
-- | the quaternion exponential function, that is, `exp <<< log` is
-- | approximately equal to the identity function. Note that `log <<< exp` is
-- | not equal to the identity function, because in general there will be many
-- | distinct quaternions `q` satisfying `exp q = p` for some fixed quaternion
-- | `p`, and the `log` function just picks one of them.
-- |
-- | This function agrees with the standard real logarithm. More generally, it
-- | agrees with the principal value of the complex logarithm, in that for a
-- | complex argument &mdash; that is, an argument with y and z components
-- | equal to zero &mdash; it returns a complex quaternion whose x component is
-- | between -π and π.
-- |
-- | Note that this function is not defined at zero.
log :: Quaternion Number -> Quaternion Number
log q =
  let
    a = realPart q
    v = vectorPart q
    normQ = norm q
    normV = Vec3.norm v
    k' = Math.acos (a / normQ)
  in
    fromReal (Math.log normQ) + fromVector (Vec3.scalarMul (k' / normV) v)

-- | The standard (Euclidean) norm of a quaternion. This makes the quaternions
-- | into a normed space; it is equivalent to the standard Euclidean norm on
-- | R^4. Defined as
-- |
-- |     norm (Quaternion w x y z) = Math.sqrt (w*w + x*x + y*y + z*z)
-- |
-- | For example:
-- |
-- |     norm (Quaternion 1.0 (-2.0) 3.0 (-4.0)) = 5.477225575051661
-- |
norm :: Quaternion Number -> Number
norm q = Math.sqrt (normSquare q)

-- | The square of the norm of a quaternion. This is slightly easier to compute
-- | than the actual norm, so may be useful in cases where you are worried
-- | about performance. Defined as
-- |
-- |     normSquare (Quaternion w x y z) = w*w + x*x + y*y + z*z
-- |
-- | For example:
-- |
-- |     normSquare (Quaternion 1.0 (-2.0) 3.0 (-4.0)) = 30.0
-- |
normSquare :: forall a. Semiring a => Quaternion a -> a
normSquare (Quaternion w x y z) = w*w + x*x + y*y + z*z

-- | The ℓ^∞ norm (when considering the quaternions as a 4-dimensional real
-- | vector space); returns the maximum absolute value of the four components.
-- | Like `norm`, this function also makes the quaternions into a normed space,
-- | although it is a little different to one to the one given by `norm`. For
-- | example:
-- |
-- |     infinityNorm (Quaternion 1.0 (-2.0) 3.0 (-4.0)) = 4.0
-- |
infinityNorm :: forall a. Ord a => Ring a => Quaternion a -> a
infinityNorm = unwrap <<< foldMap1 (Max <<< abs)

i :: forall a. Semiring a => Quaternion a
i = Quaternion zero one zero zero

j :: forall a. Semiring a => Quaternion a
j = Quaternion zero zero one zero

k :: forall a. Semiring a => Quaternion a
k = Quaternion zero zero zero one

-- | Multiplies both the real part and the vector part by the given scalar.
scalarMul :: forall a. Semiring a => a -> Quaternion a -> Quaternion a
scalarMul k' = map (k' * _)

-- | Scales the given quaternion, returning a quaternion pointing in the same
-- | direction of unit norm; multiplying this by the original quaternion's
-- | norm will give you back the original quaternion. That is, for all
-- | quaternions `q`, we have
-- |
-- |     scalarMul (norm q) (versor q) == q
-- |
versor :: Quaternion Number -> Quaternion Number
versor q = scalarMul (1.0 / norm q) q

-- | Approximate equality of quaternions, given an epsilon value specifying the
-- | maximum amount that any of the four components is allowed to differ by.
-- | Defined as
-- |
-- |     \eps p q -> infinityNorm (p - q) < eps
-- |
approxEq :: forall a. Ord a => Ring a => a -> Quaternion a -> Quaternion a -> Boolean
approxEq eps p q =
  infinityNorm (p - q) < eps
