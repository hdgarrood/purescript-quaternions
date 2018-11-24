-- | This module provides rotations in 3D space based on quaternions. It is
-- | intended to be imported qualified, for example like this:
-- |
-- |     import Data.Quaternion.Rotation (Rotation)
-- |     import Data.Quaternion.Rotation as Rotation
module Data.Quaternion.Rotation
  ( Rotation()
  , fromQuaternion
  , toQuaternion
  , fromAxisAngle
  , toAxisAngle
  , showAxisAngle
  , fromRotationMatrix
  , toRotationMatrix
  , act
  , inverse
  , normalize
  , approxEq
  ) where

import Prelude

import Data.Quaternion (Quaternion(..), conjugateBy, vectorPart, versor)
import Data.Quaternion as Q
import Data.Quaternion.Vec3 (Vec3)
import Data.Quaternion.Vec3 as Vec3
import Math as Math
import Partial (crashWith)
import Partial.Unsafe (unsafePartial)

-- | A rotation in three-dimensional space, represented by a unit quaternion
-- | (also known as a versor).
-- |
-- | The constructor is not exported to ensure that only valid rotations can
-- | be constructed.
-- |
-- | The semigroup instance provides composition of rotations; `p <> q` gives
-- | a rotation representating the rotation `q` followed by the rotation `p`.
-- | Note that in general, this is not the same as `p` followed by `q`!
newtype Rotation = Rotation (Quaternion Number)

derive instance eqRotation :: Eq Rotation

instance semigroupRotation :: Semigroup Rotation where
  append (Rotation p) (Rotation q) =
    Rotation (p * q)

instance monoidRotation :: Monoid Rotation where
  mempty = Rotation one

instance showRotation :: Show Rotation where
  show (Rotation p) = "(Rotation " <> show p <> ")"

-- | Construct a Rotation from any Quaternion, by normalizing to a versor.
fromQuaternion :: Quaternion Number -> Rotation
fromQuaternion = Rotation <<< versor

-- | Get the underlying versor.
toQuaternion :: Rotation -> Quaternion Number
toQuaternion (Rotation p) = p

-- | Approximate equality of rotations, up to a specified tolerance. Note that,
-- | for a unit quaternion `p`, the rotations represented by `p` and `-p` are
-- | identical; this function takes this into account, so that, for example:
-- |
-- |     Rotation.approxEq epsilon
-- |        (Rotation.fromQuaternion p)
-- |        (Rotation.fromQuaternion (-p))
-- |
-- | is true for all `p :: Quaternion`. Defined as
-- |
-- |     \eps p q ->
-- |        let
-- |          p' = toQuaternion p
-- |          q' = toQuaternion q
-- |        in
-- |          (Quaternion.approxEq p' q' || Quaternion.approxEq (negate p') q')
-- |
approxEq :: Number -> Rotation -> Rotation -> Boolean
approxEq eps p q =
  let
    p' = toQuaternion p
    q' = toQuaternion q
  in
    Q.approxEq eps p' q' || Q.approxEq eps (negate p') q'

-- | Construct a `Rotation` representing the rotation by the specified angle
-- | (in radians) about the specified axis. The rotation is clockwise from the
-- | point of view of someone looking along the direction of the rotation axis.
fromAxisAngle :: { angle :: Number, axis :: Vec3 Number } -> Rotation
fromAxisAngle { angle, axis } =
  let
    halfAngle = 0.5 * angle
    a = Math.sin halfAngle
  in
    unsafePartial $ case Vec3.toArray (Vec3.normalize axis) of
      [x, y, z] ->
        Rotation (Quaternion (Math.cos halfAngle) (a * x) (a * y) (a * z))

-- | Gives the angle and axis that a rotation represents. The axis returned is
-- | a unit-length vector. This is approximately an inverse of `fromAxisAngle`,
-- | in that `fromAxisAngle <<< toAxisAngle == identity`. However,
-- | `toAxisAngle <<< fromAxisAngle` is not equal to the identity function,
-- | because the axis returned will always be of unit length.
toAxisAngle :: Rotation -> { angle :: Number, axis :: Vec3 Number }
toAxisAngle (Rotation (Quaternion a b c d)) =
  let
    halfAngle = Math.acos a
    angle = halfAngle * 2.0
    k = 1.0 / Math.sin halfAngle
    axis = Vec3.vec3 (k * b) (k * c) (k * d)
  in
    { angle, axis }

-- | An alternative string representation, which can be useful for debugging.
-- | For example:
-- |
-- |     > showAxisAngle (fromQuaternion (i+j))
-- |     "Rotation.fromAxisAngle { angle: 3.141592653589793, axis: [0.7071067811865475,0.7071067811865475,0.0]}"
-- |
showAxisAngle :: Rotation -> String
showAxisAngle q =
  "Rotation.fromAxisAngle " <> show (toAxisAngle q)

-- | The inverse of a rotation; `inverse p` undoes the rotation represented by
-- | `p`. The following should hold for any rotation `p`:
-- |
-- | * `inverse p <> p == mempty`
-- | * `p <> inverse p == mempty`
-- |
-- | although note that these may only hold approximately due to floating point
-- | errors.
inverse :: Rotation -> Rotation
inverse (Rotation p) = Rotation (recip p)

-- | The action of a rotation on a vector in 3D space. This is a group action,
-- | which means that the following hold (approximately):
-- |
-- | * Identity: `act mempty == identity`
-- | * Compatibility: `act p (act q v) = act (p <> q) v`
act :: Rotation -> Vec3 Number -> Vec3 Number
act (Rotation p) v =
  unsafePartial $ case Vec3.toArray v of
    [x, y, z] ->
      vectorPart (Quaternion zero x y z `conjugateBy` p)

-- | Though all functions in this library which create a `Rotation` ensure that
-- | the underlying `Quaternion` has magnitude 1, after a sufficient number of
-- | arithmetic operations the magnitude may drift away from 1. In this case
-- | `normalize` can be used; `normalize` takes a possibly-drifted `Rotation`
-- | and rescales if it necessary, so that its magnitude returns to 1.
normalize :: Rotation -> Rotation
normalize (Rotation q) = Rotation (versor q)

-- The functions for converting to and from matrices are taken from
-- https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_.E2.86.94_quaternion

-- | Represent a Rotation as a 3-by-3 rotation matrix. The return value is an
-- | array with exactly 9 elements, in column-major order.
toRotationMatrix :: Rotation -> Array Number
toRotationMatrix (Rotation (Quaternion w x y z)) =
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
    [1.0-2.0*(yy+zz),     2.0*(xy+zw),     2.0*(xz-yw),
         2.0*(xy-zw), 1.0-2.0*(xx+zz),     2.0*(yz+xw),
         2.0*(xz+yw),     2.0*(yz-xw), 1.0-2.0*(xx+yy)]

-- | Convert a 3-by-3 rotation matrix to a Rotation representing the same
-- | rotation. The argument should be an array with exactly 9 elements, with
-- | the entries in column-major order. If the argument does not have 9
-- | elements, or if it does not represent a rotation matrix, the behaviour of
-- | this function is not defined.
fromRotationMatrix :: Partial => Array Number -> Rotation
fromRotationMatrix [a11, a21, a31, a12, a22, a32, a13, a23, a33] =
  let
    qr = 0.5 * Math.sqrt (1.0 + a11 + a22 + a33)
    qi = 0.25 / qr * (a32 - a23)
    qj = 0.25 / qr * (a13 - a31)
    qk = 0.25 / qr * (a21 - a12)
  in
    Rotation (Quaternion qr qi qj qk)
fromRotationMatrix _ =
  crashWith "Argument array is the wrong size; expected an array with 9 elements"
