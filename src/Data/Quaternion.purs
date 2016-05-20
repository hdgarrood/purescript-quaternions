module Data.Quaternion where

import Prelude
import Data.Vector3 (Vec3(), vec3)
import Math as Math

-- | A Quaternion. The type parameter denotes the underlying type. Note that
-- | the underlying type should be a reasonable approximation of the real
-- | numbers; if this is not the case, some of the functions may exhibit
-- | strange behaviour.
-- |
-- | Note that quaternion multiplication is not commutative; that is, p * q
-- | is usually not the same as q * p. Because of this, there is no Num
-- | Quaternion instance.
-- |
-- | The ModuloSemiring instance implements left division; if you need to be
-- | explicit, both types of division are also provided as separate named
-- | functions.
data Quaternion a = Quaternion a a a a

instance eqQuaternion :: Eq a => Eq (Quaternion a) where
  eq (Quaternion a1 b1 c1 d1) (Quaternion a2 b2 c2 d2) =
    a1 == a2 && b1 == b2 && c1 == c2 && d1 == d2

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

instance moduloSemiringQuaternion :: DivisionRing a => ModuloSemiring (Quaternion a) where
  div = leftDiv
  mod _ _ = zero

instance divisionRingQuaternion :: DivisionRing a => DivisionRing (Quaternion a)

realPart :: forall a. Quaternion a -> a
realPart (Quaternion a _ _ _) = a

vectorPart :: forall a. Quaternion a -> Vec3 a
vectorPart (Quaternion _ x y z) = vec3 x y z

-- | The conjugate of a quaternion. This operation negates the vector part of
-- | the quaternion.
conjugate :: forall a. Ring a => Quaternion a -> Quaternion a
conjugate (Quaternion a b c d) =
  Quaternion a (-b) (-c) (-d)

-- | The conjugate of a quaternion by another quaternion. Defined as
-- | `conjBy p q = q * p * recip q`.
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

scalarMul :: forall a. Semiring a => a -> Quaternion a -> Quaternion a
scalarMul k' (Quaternion a b c d) =
  Quaternion (k' * a) (k' * b) (k' * c) (k' * d)

-- | Gives a unit quaternion; multiplying this by the quaternion's norm will
-- | give you back the original quaternion.
versor :: Quaternion Number -> Quaternion Number
versor q = scalarMul (1.0 / norm q) q

-- | The reciprocal of a quaternion. Multiplying a quaternion by its reciprocal
-- | yields 1.
recip :: forall a. DivisionRing a => Quaternion a -> Quaternion a
recip q = scalarMul (one / normSquare q) (conjugate q)

-- | Left division; defined as `leftDiv p q = p * recip q`.
leftDiv :: forall a. DivisionRing a => Quaternion a -> Quaternion a -> Quaternion a
leftDiv p q = p * recip q

-- | Right division; defined as `rightDiv p q = recip q * p`.
rightDiv :: forall a. DivisionRing a => Quaternion a -> Quaternion a -> Quaternion a
rightDiv p q = recip q * p

-- | Approximate equality of quaternions, given an epsilon value specifying the
-- | maximum amount that any of the four components is allowed to differ by.
approxEq :: forall a. (Ord a, Ring a) => a -> Quaternion a -> Quaternion a -> Boolean
approxEq eps (Quaternion a1 b1 c1 d1) (Quaternion a2 b2 c2 d2) =
  ok a1 a2 && ok b1 b2 && ok c1 c2 && ok d1 d2
  where
  ok x y =
    let
      z = x - y
      diff = if z > zero then z else negate z
    in
      diff < eps
