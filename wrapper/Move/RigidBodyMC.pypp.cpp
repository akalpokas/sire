// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "RigidBodyMC.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/sphere.h"

#include "SireMaths/vectorproperty.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/atomselection.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculedata.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "ensemble.h"

#include "rigidbodymc.h"

#include "uniformsampler.h"

#include <QDebug>

#include <QElapsedTimer>

#include <QTime>

#include <cmath>

#include "rigidbodymc.h"

SireMove::RigidBodyMC __copy__(const SireMove::RigidBodyMC &other){ return SireMove::RigidBodyMC(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_RigidBodyMC_class(){

    { //::SireMove::RigidBodyMC
        typedef bp::class_< SireMove::RigidBodyMC, bp::bases< SireMove::MonteCarlo, SireMove::Move, SireBase::Property > > RigidBodyMC_exposer_t;
        RigidBodyMC_exposer_t RigidBodyMC_exposer = RigidBodyMC_exposer_t( "RigidBodyMC", "This class implements a rigid body Monte Carlo move that\nmay be applied to a random molecule within a MoleculeGroup\n\nAuthor: Christopher Woods\n", bp::init< bp::optional< SireBase::PropertyMap const & > >(( bp::arg("map")=SireBase::PropertyMap() ), "Null constructor") );
        bp::scope RigidBodyMC_scope( RigidBodyMC_exposer );
        RigidBodyMC_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a move that moves molecule views from the molecule group molgroup,\nselecting views randomly, with each view having an equal chance of\nbeing chosen") );
        RigidBodyMC_exposer.def( bp::init< SireMove::Sampler const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("sampler"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a move that moves molecules returned by the sampler sampler") );
        RigidBodyMC_exposer.def( bp::init< SireMove::RigidBodyMC const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::RigidBodyMC::centerOfRotation
        
            typedef ::SireMove::GetPoint const & ( ::SireMove::RigidBodyMC::*centerOfRotation_function_type)(  ) const;
            centerOfRotation_function_type centerOfRotation_function_value( &::SireMove::RigidBodyMC::centerOfRotation );
            
            RigidBodyMC_exposer.def( 
                "centerOfRotation"
                , centerOfRotation_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the function used to get the center of rotation of each molecule" );
        
        }
        { //::SireMove::RigidBodyMC::disableReflectionSphere
        
            typedef void ( ::SireMove::RigidBodyMC::*disableReflectionSphere_function_type)(  ) ;
            disableReflectionSphere_function_type disableReflectionSphere_function_value( &::SireMove::RigidBodyMC::disableReflectionSphere );
            
            RigidBodyMC_exposer.def( 
                "disableReflectionSphere"
                , disableReflectionSphere_function_value
                , "Completely switch off use of the reflection sphere or volume" );
        
        }
        { //::SireMove::RigidBodyMC::disableReflectionVolume
        
            typedef void ( ::SireMove::RigidBodyMC::*disableReflectionVolume_function_type)(  ) ;
            disableReflectionVolume_function_type disableReflectionVolume_function_value( &::SireMove::RigidBodyMC::disableReflectionVolume );
            
            RigidBodyMC_exposer.def( 
                "disableReflectionVolume"
                , disableReflectionVolume_function_value
                , "Completely switch off use of the reflection sphere or volume" );
        
        }
        { //::SireMove::RigidBodyMC::extract
        
            typedef ::SireMol::Molecules ( ::SireMove::RigidBodyMC::*extract_function_type)( ::SireMol::Molecules const & ) const;
            extract_function_type extract_function_value( &::SireMove::RigidBodyMC::extract );
            
            RigidBodyMC_exposer.def( 
                "extract"
                , extract_function_value
                , ( bp::arg("molecules") )
                , "Extract from mols all molecules that are within the reflection sphere  volume\ndescribed by this move. This allows the molecules that will be affected by this\nmove to be separated out from the rest.\nNote that this returns all molecules if a reflection sphere or volume is not used\n" );
        
        }
        { //::SireMove::RigidBodyMC::extract
        
            typedef ::SireMol::Molecules ( ::SireMove::RigidBodyMC::*extract_function_type)( ::SireMol::Molecules const &,::SireUnits::Dimension::Length ) const;
            extract_function_type extract_function_value( &::SireMove::RigidBodyMC::extract );
            
            RigidBodyMC_exposer.def( 
                "extract"
                , extract_function_value
                , ( bp::arg("molecules"), bp::arg("buffer") )
                , "Extract from mols all molecules that are within the reflection sphere  volume\ndescribed by this move. This allows the molecules that will be affected by this\nmove to be separated out from the rest.\nIf buffer is greater than zero, then this will also extract molecules that are\nwithin buffer of the reflection volume, and will additionally translate those\nmolecules so that they are moved into the volume. The buffer can be used\nto pack more molecules into the volume than would be achieved just with a simple\nselection.\nNote that this returns all molecules if a reflection sphere or volume is not used\n" );
        
        }
        { //::SireMove::RigidBodyMC::maximumRotation
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMove::RigidBodyMC::*maximumRotation_function_type)(  ) const;
            maximumRotation_function_type maximumRotation_function_value( &::SireMove::RigidBodyMC::maximumRotation );
            
            RigidBodyMC_exposer.def( 
                "maximumRotation"
                , maximumRotation_function_value
                , "Return the maximum rotation for each move" );
        
        }
        { //::SireMove::RigidBodyMC::maximumTranslation
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::RigidBodyMC::*maximumTranslation_function_type)(  ) const;
            maximumTranslation_function_type maximumTranslation_function_value( &::SireMove::RigidBodyMC::maximumTranslation );
            
            RigidBodyMC_exposer.def( 
                "maximumTranslation"
                , maximumTranslation_function_value
                , "Return the maximum translation for each move" );
        
        }
        { //::SireMove::RigidBodyMC::moleculeGroup
        
            typedef ::SireMol::MoleculeGroup const & ( ::SireMove::RigidBodyMC::*moleculeGroup_function_type)(  ) const;
            moleculeGroup_function_type moleculeGroup_function_value( &::SireMove::RigidBodyMC::moleculeGroup );
            
            RigidBodyMC_exposer.def( 
                "moleculeGroup"
                , moleculeGroup_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the molecule group from which random molecule views\nare drawn" );
        
        }
        { //::SireMove::RigidBodyMC::move
        
            typedef void ( ::SireMove::RigidBodyMC::*move_function_type)( ::SireSystem::System &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::RigidBodyMC::move );
            
            RigidBodyMC_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Attempt n rigid body moves of the views of the system system" );
        
        }
        RigidBodyMC_exposer.def( bp::self != bp::self );
        { //::SireMove::RigidBodyMC::operator=
        
            typedef ::SireMove::RigidBodyMC & ( ::SireMove::RigidBodyMC::*assign_function_type)( ::SireMove::RigidBodyMC const & ) ;
            assign_function_type assign_function_value( &::SireMove::RigidBodyMC::operator= );
            
            RigidBodyMC_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        RigidBodyMC_exposer.def( bp::self == bp::self );
        { //::SireMove::RigidBodyMC::reflectedVolume
        
            typedef double ( ::SireMove::RigidBodyMC::*reflectedVolume_function_type)(  ) const;
            reflectedVolume_function_type reflectedVolume_function_value( &::SireMove::RigidBodyMC::reflectedVolume );
            
            RigidBodyMC_exposer.def( 
                "reflectedVolume"
                , reflectedVolume_function_value
                , "Return the volume of space occupied by the reflection volume" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionSphereCenter
        
            typedef ::SireMaths::Vector ( ::SireMove::RigidBodyMC::*reflectionSphereCenter_function_type)(  ) const;
            reflectionSphereCenter_function_type reflectionSphereCenter_function_value( &::SireMove::RigidBodyMC::reflectionSphereCenter );
            
            RigidBodyMC_exposer.def( 
                "reflectionSphereCenter"
                , reflectionSphereCenter_function_value
                , "Return the center of the reflection sphere. Returns a zero vector\nif a reflection sphere is not being used" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionSphereCenter
        
            typedef ::SireMaths::Vector ( ::SireMove::RigidBodyMC::*reflectionSphereCenter_function_type)( ::SireMol::MolNum ) const;
            reflectionSphereCenter_function_type reflectionSphereCenter_function_value( &::SireMove::RigidBodyMC::reflectionSphereCenter );
            
            RigidBodyMC_exposer.def( 
                "reflectionSphereCenter"
                , reflectionSphereCenter_function_value
                , ( bp::arg("molnum") )
                , "Return the center of the reflection sphere for molecule molnum. Returns a zero vector\nif a reflection sphere is not being used" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionSphereCenter
        
            typedef ::SireMaths::Vector ( ::SireMove::RigidBodyMC::*reflectionSphereCenter_function_type)( ::SireMol::MoleculeView const & ) const;
            reflectionSphereCenter_function_type reflectionSphereCenter_function_value( &::SireMove::RigidBodyMC::reflectionSphereCenter );
            
            RigidBodyMC_exposer.def( 
                "reflectionSphereCenter"
                , reflectionSphereCenter_function_value
                , ( bp::arg("molview") )
                , "Return the center of the reflection sphere for molecule mol. Returns a null vector\nif a reflection sphere is not being used" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionSphereRadius
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::RigidBodyMC::*reflectionSphereRadius_function_type)(  ) const;
            reflectionSphereRadius_function_type reflectionSphereRadius_function_value( &::SireMove::RigidBodyMC::reflectionSphereRadius );
            
            RigidBodyMC_exposer.def( 
                "reflectionSphereRadius"
                , reflectionSphereRadius_function_value
                , "Return the radius of the reflection sphere. This returns zero\nif the reflection sphere is not being used" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionSphereRadius
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::RigidBodyMC::*reflectionSphereRadius_function_type)( ::SireMol::MolNum ) const;
            reflectionSphereRadius_function_type reflectionSphereRadius_function_value( &::SireMove::RigidBodyMC::reflectionSphereRadius );
            
            RigidBodyMC_exposer.def( 
                "reflectionSphereRadius"
                , reflectionSphereRadius_function_value
                , ( bp::arg("molnum") )
                , "Return the radius of the reflection sphere for molecule molnum. This returns zero\nif the reflection sphere is not being used" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionSphereRadius
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::RigidBodyMC::*reflectionSphereRadius_function_type)( ::SireMol::MoleculeView const & ) const;
            reflectionSphereRadius_function_type reflectionSphereRadius_function_value( &::SireMove::RigidBodyMC::reflectionSphereRadius );
            
            RigidBodyMC_exposer.def( 
                "reflectionSphereRadius"
                , reflectionSphereRadius_function_value
                , ( bp::arg("molview") )
                , "Return the radius of the reflection sphere for molecule mol. This returns zero\nif the reflection sphere is not being used" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionVolume
        
            typedef ::QVector< SireMaths::Sphere > ( ::SireMove::RigidBodyMC::*reflectionVolume_function_type)(  ) const;
            reflectionVolume_function_type reflectionVolume_function_value( &::SireMove::RigidBodyMC::reflectionVolume );
            
            RigidBodyMC_exposer.def( 
                "reflectionVolume"
                , reflectionVolume_function_value
                , "Return the set of spheres that are used to define the reflection volume.\nThis returns an empty list if the reflection volume  sphere is not used" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionVolumePoints
        
            typedef ::QVector< SireMaths::Vector > ( ::SireMove::RigidBodyMC::*reflectionVolumePoints_function_type)(  ) const;
            reflectionVolumePoints_function_type reflectionVolumePoints_function_value( &::SireMove::RigidBodyMC::reflectionVolumePoints );
            
            RigidBodyMC_exposer.def( 
                "reflectionVolumePoints"
                , reflectionVolumePoints_function_value
                , "Return the points used to define the reflection volume. This will\nbe a single point if a reflection sphere is used. If a reflection\nsphere or volume is not used then this returns an empty list" );
        
        }
        { //::SireMove::RigidBodyMC::reflectionVolumeRadius
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::RigidBodyMC::*reflectionVolumeRadius_function_type)(  ) const;
            reflectionVolumeRadius_function_type reflectionVolumeRadius_function_value( &::SireMove::RigidBodyMC::reflectionVolumeRadius );
            
            RigidBodyMC_exposer.def( 
                "reflectionVolumeRadius"
                , reflectionVolumeRadius_function_value
                , "Return the reflection volume radius (same as the reflection sphere radius)" );
        
        }
        { //::SireMove::RigidBodyMC::sampler
        
            typedef ::SireMove::Sampler const & ( ::SireMove::RigidBodyMC::*sampler_function_type)(  ) const;
            sampler_function_type sampler_function_value( &::SireMove::RigidBodyMC::sampler );
            
            RigidBodyMC_exposer.def( 
                "sampler"
                , sampler_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the sampler that is used to draw random molecules" );
        
        }
        { //::SireMove::RigidBodyMC::setCenterOfRotation
        
            typedef void ( ::SireMove::RigidBodyMC::*setCenterOfRotation_function_type)( ::SireMove::GetPoint const & ) ;
            setCenterOfRotation_function_type setCenterOfRotation_function_value( &::SireMove::RigidBodyMC::setCenterOfRotation );
            
            RigidBodyMC_exposer.def( 
                "setCenterOfRotation"
                , setCenterOfRotation_function_value
                , ( bp::arg("center_function") )
                , "Set the function used to get the center of rotation for each molecule" );
        
        }
        { //::SireMove::RigidBodyMC::setGenerator
        
            typedef void ( ::SireMove::RigidBodyMC::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::RigidBodyMC::setGenerator );
            
            RigidBodyMC_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") )
                , "Set the random number generator used by this move" );
        
        }
        { //::SireMove::RigidBodyMC::setMaximumRotation
        
            typedef void ( ::SireMove::RigidBodyMC::*setMaximumRotation_function_type)( ::SireUnits::Dimension::Angle ) ;
            setMaximumRotation_function_type setMaximumRotation_function_value( &::SireMove::RigidBodyMC::setMaximumRotation );
            
            RigidBodyMC_exposer.def( 
                "setMaximumRotation"
                , setMaximumRotation_function_value
                , ( bp::arg("max_rotation") )
                , "Set the maximum delta for any rotation" );
        
        }
        { //::SireMove::RigidBodyMC::setMaximumTranslation
        
            typedef void ( ::SireMove::RigidBodyMC::*setMaximumTranslation_function_type)( ::SireUnits::Dimension::Length ) ;
            setMaximumTranslation_function_type setMaximumTranslation_function_value( &::SireMove::RigidBodyMC::setMaximumTranslation );
            
            RigidBodyMC_exposer.def( 
                "setMaximumTranslation"
                , setMaximumTranslation_function_value
                , ( bp::arg("max_translation") )
                , "Set the maximum delta for any translation" );
        
        }
        { //::SireMove::RigidBodyMC::setReflectionSphere
        
            typedef void ( ::SireMove::RigidBodyMC::*setReflectionSphere_function_type)( ::SireMaths::Vector,::SireUnits::Dimension::Length ) ;
            setReflectionSphere_function_type setReflectionSphere_function_value( &::SireMove::RigidBodyMC::setReflectionSphere );
            
            RigidBodyMC_exposer.def( 
                "setReflectionSphere"
                , setReflectionSphere_function_value
                , ( bp::arg("sphere_center"), bp::arg("sphere_radius") )
                , "Turn on rigid body move reflections. This makes sure that only\nmolecules within the specified sphere can be moved, and any moves\nare reflected so that these molecules will always remain within\nthe sphere" );
        
        }
        { //::SireMove::RigidBodyMC::setReflectionSphere
        
            typedef void ( ::SireMove::RigidBodyMC::*setReflectionSphere_function_type)( ::SireMol::MolNum,::SireMaths::Vector,::SireUnits::Dimension::Length ) ;
            setReflectionSphere_function_type setReflectionSphere_function_value( &::SireMove::RigidBodyMC::setReflectionSphere );
            
            RigidBodyMC_exposer.def( 
                "setReflectionSphere"
                , setReflectionSphere_function_value
                , ( bp::arg("molnum"), bp::arg("sphere_center"), bp::arg("sphere_radius") )
                , "Turn on rigid body move reflections for molecule molnum. This makes sure\nthis molecule is moved only within the sphere centered at sphere_center\nwith radius sphere_radius. Any moves are reflected so that\nthis molecule will always remain within the sphere" );
        
        }
        { //::SireMove::RigidBodyMC::setReflectionSphere
        
            typedef void ( ::SireMove::RigidBodyMC::*setReflectionSphere_function_type)( ::SireMol::MoleculeView const &,::SireMaths::Vector,::SireUnits::Dimension::Length ) ;
            setReflectionSphere_function_type setReflectionSphere_function_value( &::SireMove::RigidBodyMC::setReflectionSphere );
            
            RigidBodyMC_exposer.def( 
                "setReflectionSphere"
                , setReflectionSphere_function_value
                , ( bp::arg("molview"), bp::arg("sphere_center"), bp::arg("sphere_radius") )
                , "Turn on rigid body move reflections for molecule mol. This makes sure\nthis molecule is moved only within the sphere centered at sphere_center\nwith radius sphere_radius. Any moves are reflected so that\nthis molecule will always remain within the sphere" );
        
        }
        { //::SireMove::RigidBodyMC::setReflectionVolume
        
            typedef void ( ::SireMove::RigidBodyMC::*setReflectionVolume_function_type)( ::QVector< SireMaths::Vector > const &,::SireUnits::Dimension::Length ) ;
            setReflectionVolume_function_type setReflectionVolume_function_value( &::SireMove::RigidBodyMC::setReflectionVolume );
            
            RigidBodyMC_exposer.def( 
                "setReflectionVolume"
                , setReflectionVolume_function_value
                , ( bp::arg("points"), bp::arg("radius") )
                , "Turn on and specify a reflection volume. This is like a reflection sphere,\nexcept now the reflection volume is formed as the intersection of the spheres\nwhose centers are in points, all with radii radius. This replaces any\nreflection sphere set (although not molecule-specific reflection spheres)" );
        
        }
        { //::SireMove::RigidBodyMC::setReflectionVolume
        
            typedef void ( ::SireMove::RigidBodyMC::*setReflectionVolume_function_type)( ::SireMol::MoleculeView const &,::SireUnits::Dimension::Length,bool,::SireBase::PropertyMap const & ) ;
            setReflectionVolume_function_type setReflectionVolume_function_value( &::SireMove::RigidBodyMC::setReflectionVolume );
            
            RigidBodyMC_exposer.def( 
                "setReflectionVolume"
                , setReflectionVolume_function_value
                , ( bp::arg("molecule"), bp::arg("radius"), bp::arg("heavy_atoms_only"), bp::arg("map")=SireBase::PropertyMap() )
                , "Turn on and specify the reflection volume by using the coordinates of\nall of the atoms in the passed molecule view (excluding light atoms\nif heavy_atoms_only is true)" );
        
        }
        { //::SireMove::RigidBodyMC::setReflectionVolume
        
            typedef void ( ::SireMove::RigidBodyMC::*setReflectionVolume_function_type)( ::SireMol::MoleculeView const &,::SireUnits::Dimension::Length,::SireBase::PropertyMap const & ) ;
            setReflectionVolume_function_type setReflectionVolume_function_value( &::SireMove::RigidBodyMC::setReflectionVolume );
            
            RigidBodyMC_exposer.def( 
                "setReflectionVolume"
                , setReflectionVolume_function_value
                , ( bp::arg("molecule"), bp::arg("radius"), bp::arg("map")=SireBase::PropertyMap() )
                , "Turn on and specify the reflection volume by using the coordinates of\nall of the heavy atoms in the passed molecule view" );
        
        }
        { //::SireMove::RigidBodyMC::setSampler
        
            typedef void ( ::SireMove::RigidBodyMC::*setSampler_function_type)( ::SireMove::Sampler const & ) ;
            setSampler_function_type setSampler_function_value( &::SireMove::RigidBodyMC::setSampler );
            
            RigidBodyMC_exposer.def( 
                "setSampler"
                , setSampler_function_value
                , ( bp::arg("sampler") )
                , "Set the sampler (and contained molecule group) that provides\nthe random molecules to be moved. This gives the sampler the\nsame random number generator that is used by this move" );
        
        }
        { //::SireMove::RigidBodyMC::setSampler
        
            typedef void ( ::SireMove::RigidBodyMC::*setSampler_function_type)( ::SireMol::MoleculeGroup const & ) ;
            setSampler_function_type setSampler_function_value( &::SireMove::RigidBodyMC::setSampler );
            
            RigidBodyMC_exposer.def( 
                "setSampler"
                , setSampler_function_value
                , ( bp::arg("molgroup") )
                , "Set the sampler to be one that selects views at random\nfrom the molecule group molgroup (each view has an\nequal chance of being chosen)" );
        
        }
        { //::SireMove::RigidBodyMC::setSharedRotationCenter
        
            typedef void ( ::SireMove::RigidBodyMC::*setSharedRotationCenter_function_type)( bool ) ;
            setSharedRotationCenter_function_type setSharedRotationCenter_function_value( &::SireMove::RigidBodyMC::setSharedRotationCenter );
            
            RigidBodyMC_exposer.def( 
                "setSharedRotationCenter"
                , setSharedRotationCenter_function_value
                , ( bp::arg("on") )
                , "Set whether or not to use the same rotation center for all\nsynchronised molecules" );
        
        }
        { //::SireMove::RigidBodyMC::setSynchronisedRotation
        
            typedef void ( ::SireMove::RigidBodyMC::*setSynchronisedRotation_function_type)( bool ) ;
            setSynchronisedRotation_function_type setSynchronisedRotation_function_value( &::SireMove::RigidBodyMC::setSynchronisedRotation );
            
            RigidBodyMC_exposer.def( 
                "setSynchronisedRotation"
                , setSynchronisedRotation_function_value
                , ( bp::arg("on") )
                , "Set whether or not to synchronise rotation of all of the views" );
        
        }
        { //::SireMove::RigidBodyMC::setSynchronisedTranslation
        
            typedef void ( ::SireMove::RigidBodyMC::*setSynchronisedTranslation_function_type)( bool ) ;
            setSynchronisedTranslation_function_type setSynchronisedTranslation_function_value( &::SireMove::RigidBodyMC::setSynchronisedTranslation );
            
            RigidBodyMC_exposer.def( 
                "setSynchronisedTranslation"
                , setSynchronisedTranslation_function_value
                , ( bp::arg("on") )
                , "Set whether or not to synchronise translation of all of the views" );
        
        }
        { //::SireMove::RigidBodyMC::sharedRotationCenter
        
            typedef bool ( ::SireMove::RigidBodyMC::*sharedRotationCenter_function_type)(  ) const;
            sharedRotationCenter_function_type sharedRotationCenter_function_value( &::SireMove::RigidBodyMC::sharedRotationCenter );
            
            RigidBodyMC_exposer.def( 
                "sharedRotationCenter"
                , sharedRotationCenter_function_value
                , "Return whether or not synchronised rotation uses the same\ncenter of rotation for all molecules" );
        
        }
        { //::SireMove::RigidBodyMC::synchronisedRotation
        
            typedef bool ( ::SireMove::RigidBodyMC::*synchronisedRotation_function_type)(  ) const;
            synchronisedRotation_function_type synchronisedRotation_function_value( &::SireMove::RigidBodyMC::synchronisedRotation );
            
            RigidBodyMC_exposer.def( 
                "synchronisedRotation"
                , synchronisedRotation_function_value
                , "Return whether or not rotation of all molecules is synchronised" );
        
        }
        { //::SireMove::RigidBodyMC::synchronisedTranslation
        
            typedef bool ( ::SireMove::RigidBodyMC::*synchronisedTranslation_function_type)(  ) const;
            synchronisedTranslation_function_type synchronisedTranslation_function_value( &::SireMove::RigidBodyMC::synchronisedTranslation );
            
            RigidBodyMC_exposer.def( 
                "synchronisedTranslation"
                , synchronisedTranslation_function_value
                , "Return whether or not translation of all molecules is synchronised" );
        
        }
        { //::SireMove::RigidBodyMC::toString
        
            typedef ::QString ( ::SireMove::RigidBodyMC::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::RigidBodyMC::toString );
            
            RigidBodyMC_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this move" );
        
        }
        { //::SireMove::RigidBodyMC::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::RigidBodyMC::typeName );
            
            RigidBodyMC_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMove::RigidBodyMC::usesReflectionMoves
        
            typedef bool ( ::SireMove::RigidBodyMC::*usesReflectionMoves_function_type)(  ) const;
            usesReflectionMoves_function_type usesReflectionMoves_function_value( &::SireMove::RigidBodyMC::usesReflectionMoves );
            
            RigidBodyMC_exposer.def( 
                "usesReflectionMoves"
                , usesReflectionMoves_function_value
                , "Return whether or not these moves use a reflection sphere" );
        
        }
        { //::SireMove::RigidBodyMC::usesReflectionMoves
        
            typedef bool ( ::SireMove::RigidBodyMC::*usesReflectionMoves_function_type)( ::SireMol::MolNum ) const;
            usesReflectionMoves_function_type usesReflectionMoves_function_value( &::SireMove::RigidBodyMC::usesReflectionMoves );
            
            RigidBodyMC_exposer.def( 
                "usesReflectionMoves"
                , usesReflectionMoves_function_value
                , ( bp::arg("molnum") )
                , "Return whether or not these moves use a reflection sphere for molecule molnum" );
        
        }
        { //::SireMove::RigidBodyMC::usesReflectionMoves
        
            typedef bool ( ::SireMove::RigidBodyMC::*usesReflectionMoves_function_type)( ::SireMol::MoleculeView const & ) const;
            usesReflectionMoves_function_type usesReflectionMoves_function_value( &::SireMove::RigidBodyMC::usesReflectionMoves );
            
            RigidBodyMC_exposer.def( 
                "usesReflectionMoves"
                , usesReflectionMoves_function_value
                , ( bp::arg("molview") )
                , "Return whether or not these moves use a reflection sphere for molecule mol" );
        
        }
        { //::SireMove::RigidBodyMC::usesReflectionVolume
        
            typedef bool ( ::SireMove::RigidBodyMC::*usesReflectionVolume_function_type)(  ) const;
            usesReflectionVolume_function_type usesReflectionVolume_function_value( &::SireMove::RigidBodyMC::usesReflectionVolume );
            
            RigidBodyMC_exposer.def( 
                "usesReflectionVolume"
                , usesReflectionVolume_function_value
                , "Return whether or not rigid body moves are confined to a reflection volume.\nNote that this also returns true if we are using a reflection sphere\n(as a sphere is just a special case of a reflection volume)" );
        
        }
        RigidBodyMC_exposer.staticmethod( "typeName" );
        RigidBodyMC_exposer.def( "__copy__", &__copy__);
        RigidBodyMC_exposer.def( "__deepcopy__", &__copy__);
        RigidBodyMC_exposer.def( "clone", &__copy__);
        RigidBodyMC_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::RigidBodyMC >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RigidBodyMC_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::RigidBodyMC >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RigidBodyMC_exposer.def( "__str__", &__str__< ::SireMove::RigidBodyMC > );
        RigidBodyMC_exposer.def( "__repr__", &__str__< ::SireMove::RigidBodyMC > );
    }

}
