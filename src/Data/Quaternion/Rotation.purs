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
  , normalize
  , toMat4
  ) where

import Prelude
import Data.Monoid (class Monoid)
import Data.Vector as Vec
import Data.Vector3 (Vec3)
import Data.Vector3 as Vec3
import Data.Matrix4 as Mat4
import Data.Matrix4 (Mat4)

import Math as Math
import Partial.Unsafe (unsafePartial)

import Data.Quaternion (Quaternion(..), conjugateBy, vectorPart, versor)

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
    unsafePartial $ case Vec.toArray (Vec.normalize axis) of
      [x, y, z] ->
        Rotation (Quaternion (Math.cos halfAngle) (a * x) (a * y) (a * z))

-- | Gives the angle and axis that a rotation represents. This is the inverse
-- | of `fromAngleAxis`.
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

-- | An alternative string representation, which can be useful for debugging.
showAngleAxis :: Rotation Number -> String
showAngleAxis q =
  case toAngleAxis q of
    { angle, axis } ->
      "Rotation.fromAngleAxis " <>
       "{ angle: " <> show angle <>
       ", axis: " <> show axis <> "}"

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
  unsafePartial $ case Vec.toArray v of
    [x, y, z] ->
      vectorPart (Quaternion zero x y z `conjugateBy` p)

-- | Though all functions in the library that create a Rotation ensure that
-- | the underlying Quaternion has unit magnitude, some computations may result
-- | in the magnitude drifting away from 1.0 for numerical or other reasons.
-- | normalize takes a possibly-drifted Rotation and returns a proper Rotation.
normalize :: Rotation Number -> Rotation Number
normalize (Rotation q) = Rotation (versor q)

toMat4 :: Rotation Number -> Mat4
toMat4 (Rotation (Quaternion w x y z) ) =
  let
    xx = x * x
    xy = x * y
    xz = x * z
    xw = x * w
    yy = y * y
    yz = y * z
    yw = y * w
    zz = z * z
    zw = z * w
  in
    Mat4.mat4 [1.0-2.0*(yy+zz),     2.0*(xy+zw),     2.0*(xz-yw), 0.0,
                   2.0*(xy-zw), 1.0-2.0*(xx+zz),     2.0*(yz+xw), 0.0,
                   2.0*(xz+yw),     2.0*(yz-xw), 1.0-2.0*(xx+yy), 0.0,
                           0.0,             0.0,             0.0, 1.0]
