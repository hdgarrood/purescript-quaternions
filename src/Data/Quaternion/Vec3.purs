-- | A minimal type for vectors in 3D space. This type is a newtype over `Array
-- | a`, with constructors to ensure that we do in fact have exactly three
-- | elements.
module Data.Quaternion.Vec3 (Vec3, vec3, toArray, fromArray, normalize) where

import Prelude

import Math as Math
import Partial.Unsafe (unsafeCrashWith)

newtype Vec3 a = Vec3 (Array a)

derive instance eqVec3 :: Eq a => Eq (Vec3 a)
derive instance ordVec3 :: Ord a => Ord (Vec3 a)

instance showVec3 :: Show a => Show (Vec3 a) where
  show (Vec3 v) = show v

-- | Construct a `Vec3` from three arguments.
vec3 :: forall a. a -> a -> a -> Vec3 a
vec3 x y z = Vec3 [x, y, z]

-- | Convert a `Vec3` to a three-element array.
toArray :: forall a. Vec3 a -> Array a
toArray (Vec3 v) = v

-- | Convert a three-element array to a `Vec3`. If the argument does not have
-- | three elements, the behaviour of this function is undefined.
fromArray :: forall a. Partial => Array a -> Vec3 a
fromArray = Vec3

-- | Normalize a vector, returning a unit vector pointing in the same
-- | direction. Attempting to normalize the zero vector simply returns the zero
-- | vector.
normalize :: Vec3 Number -> Vec3 Number
normalize (Vec3 [x,y,z]) =
  let
    k = Math.sqrt (x*x + y*y + z*z)
  in
    Vec3 [k*x, k*y, k*z]
normalize _ =
  unsafeCrashWith "Vec3 invariant violated"
