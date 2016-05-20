-- | This module provides rotations in 3D space based on quaternions. It is
-- | intended to be imported qualified, for example like this:
-- |
-- |     import Data.Quaternion.Rotation (Rotation)
-- |     import Data.Quaternion.Rotation as Rotation
module Data.Quaternion.Rotation
  ( Rotation()
  , fromQuaternion
  , toQuaternion
  , fromAngleAxis
  , toAngleAxis
  , showAngleAxis
  , inverse
  , act
  ) where

import Prelude
import Data.Monoid (class Monoid)
import Data.Vector as Vec
import Data.Vector3 (Vec3)
import Data.Vector3 as Vec3
import Math as Math

import Data.Quaternion (Quaternion(..), conjugateBy, vectorPart, recip, versor)

-- | A rotation in three-dimensional space, represented by a unit quaternion
-- | (also known as a versor).
-- |
-- | The constructor is not exported to ensure that only valid rotations can
-- | be constructed.
-- |
-- | The semigroup instance provides composition of rotations; `p <> q` gives
-- | a rotation representating the rotation `q` followed by the rotation `p`.
-- | Note that in general, this is not the same as `p` followed by `q`!
newtype Rotation a = Rotation (Quaternion a)

-- | Construct a Rotation from any Quaternion, by normalizing to a versor.
fromQuaternion :: Quaternion Number -> Rotation Number
fromQuaternion = Rotation <<< versor

-- | Get the underlying versor.
toQuaternion :: forall a. Rotation a -> Quaternion a
toQuaternion (Rotation p) = p

-- | Construct a `Rotation` representing the rotation by the specified angle
-- | (in radians) about the specified axis. The rotation is clockwise from the
-- | point of view of someone looking along the direction of the rotation axis.
fromAngleAxis :: { angle :: Number, axis :: Vec3 Number } -> Rotation Number
fromAngleAxis { angle, axis } =
  let
    halfAngle = 0.5 * angle
    a = Math.sin halfAngle
  in
    case Vec.toArray (Vec.normalize axis) of
      [x, y, z] ->
        Rotation (Quaternion (Math.cos halfAngle) (a * x) (a * y) (a * z))

-- | Gives the angle and axis that a rotation represents. This is the inverse
-- | of `construct`.
toAngleAxis :: Rotation Number -> { angle :: Number, axis :: Vec3 Number }
toAngleAxis (Rotation (Quaternion a b c d)) =
  let
    halfAngle = Math.acos a
    angle = halfAngle * 2.0
    k = 1.0 / Math.sin halfAngle
    axis = Vec3.vec3 (k * b) (k * c) (k * d)
  in
    { angle, axis }

instance semigroupRotation :: Ring a => Semigroup (Rotation a) where
  append (Rotation p) (Rotation q) =
    Rotation (p * q)

instance monoidRotation :: Ring a => Monoid (Rotation a) where
  mempty = Rotation one

instance eqRotation :: Eq a => Eq (Rotation a) where
  eq (Rotation p) (Rotation q) = p == q

instance showRotation :: Show a => Show (Rotation a) where
  show (Rotation p) = "(Rotation " <> show p <> ")"

-- | An alternative string representation for debugging.
showAngleAxis :: Rotation Number -> String
showAngleAxis q =
  case toAngleAxis q of
    { angle, axis } ->
      "(Rotation: angle=" <> show angle <> " axis=" <> show (Vec.toArray axis) <> ")"

-- | The inverse of a rotation. The following should hold for any rotation `p`:
-- |
-- | * `inverse p <> p == mempty`
-- | * `p <> inverse p == mempty`
inverse :: forall a. DivisionRing a => Rotation a -> Rotation a
inverse (Rotation p) = Rotation (recip p)

-- | The action of a rotation on a vector in 3D space. This is a group action,
-- | which means that the following hold:
-- |
-- | * Identity: `act mempty == id`
-- | * Compatibility: `act p (act q v) = act (p <> q) v`
act :: forall a. DivisionRing a => Rotation a -> Vec3 a -> Vec3 a
act (Rotation p) v =
  case Vec.toArray v of
    [x, y, z] ->
      vectorPart (Quaternion zero x y z `conjugateBy` p)
