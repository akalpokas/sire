// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "BoxPatching.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "patching.h"

#include "patching.h"

SireVol::BoxPatching __copy__(const SireVol::BoxPatching &other){ return SireVol::BoxPatching(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_BoxPatching_class(){

    { //::SireVol::BoxPatching
        typedef bp::class_< SireVol::BoxPatching, bp::bases< SireVol::Patching, SireBase::Property > > BoxPatching_exposer_t;
        BoxPatching_exposer_t BoxPatching_exposer = BoxPatching_exposer_t( "BoxPatching", "This is a simple patching scheme that divides space up into a\nseries of cuboidal boxes (3D grid)\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope BoxPatching_scope( BoxPatching_exposer );
        BoxPatching_exposer.def( bp::init< SireVol::Space const & >(( bp::arg("space") ), "Construct for the passed space - this tries to divide the space\nusing a patch size of 8 A, placing the center of the grid at (0,0,0)") );
        BoxPatching_exposer.def( bp::init< SireVol::Space const &, SireMaths::Vector const & >(( bp::arg("space"), bp::arg("center") ), "Construct for the passed space, placing the center of the grid at center.\nThis tries to divide space using a patch size of 8 A") );
        BoxPatching_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length >(( bp::arg("space"), bp::arg("patch_size") ), "Construct for the passed space, using the passed patch size. This will try\nto build a cubic grid of patches where the grid dimension is approximately\npatch_size, with the center of the grid at (0,0,0)") );
        BoxPatching_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length, SireMaths::Vector const & >(( bp::arg("space"), bp::arg("patch_size"), bp::arg("center") ), "Construct for the passed space, placing the center of the grid at center,\ntrying to construct a grid that divides space using a patch size of\napproximately patch_size\nNote that this patching is only compatible with cartesian spaces.\nThrow: SireError::incompatible_error\n") );
        BoxPatching_exposer.def( bp::init< SireVol::BoxPatching const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireVol::BoxPatching::center
        
            typedef ::SireMaths::Vector ( ::SireVol::BoxPatching::*center_function_type)(  ) const;
            center_function_type center_function_value( &::SireVol::BoxPatching::center );
            
            BoxPatching_exposer.def( 
                "center"
                , center_function_value
                , "Return the center of the patching grid" );
        
        }
        { //::SireVol::BoxPatching::nPatches
        
            typedef int ( ::SireVol::BoxPatching::*nPatches_function_type)(  ) const;
            nPatches_function_type nPatches_function_value( &::SireVol::BoxPatching::nPatches );
            
            BoxPatching_exposer.def( 
                "nPatches"
                , nPatches_function_value
                , "Return the number of patches" );
        
        }
        BoxPatching_exposer.def( bp::self != bp::self );
        { //::SireVol::BoxPatching::operator=
        
            typedef ::SireVol::BoxPatching & ( ::SireVol::BoxPatching::*assign_function_type)( ::SireVol::BoxPatching const & ) ;
            assign_function_type assign_function_value( &::SireVol::BoxPatching::operator= );
            
            BoxPatching_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        BoxPatching_exposer.def( bp::self == bp::self );
        { //::SireVol::BoxPatching::patchBox
        
            typedef ::SireVol::AABox ( ::SireVol::BoxPatching::*patchBox_function_type)( int ) const;
            patchBox_function_type patchBox_function_value( &::SireVol::BoxPatching::patchBox );
            
            BoxPatching_exposer.def( 
                "patchBox"
                , patchBox_function_value
                , ( bp::arg("i") )
                , "Return the AABox that completely encloses the ith patch box" );
        
        }
        { //::SireVol::BoxPatching::patchBox
        
            typedef ::SireVol::AABox ( ::SireVol::BoxPatching::*patchBox_function_type)( ::SireMaths::Vector const & ) const;
            patchBox_function_type patchBox_function_value( &::SireVol::BoxPatching::patchBox );
            
            BoxPatching_exposer.def( 
                "patchBox"
                , patchBox_function_value
                , ( bp::arg("point") )
                , "Return the AABox of that completely encloses the patch that contains the\npoint point" );
        
        }
        { //::SireVol::BoxPatching::patchDimension
        
            typedef ::SireMaths::Vector ( ::SireVol::BoxPatching::*patchDimension_function_type)(  ) const;
            patchDimension_function_type patchDimension_function_value( &::SireVol::BoxPatching::patchDimension );
            
            BoxPatching_exposer.def( 
                "patchDimension"
                , patchDimension_function_value
                , "Return the dimensions of each path (the lengths of each side of the box)" );
        
        }
        { //::SireVol::BoxPatching::patchIndex
        
            typedef int ( ::SireVol::BoxPatching::*patchIndex_function_type)( ::SireMaths::Vector const & ) const;
            patchIndex_function_type patchIndex_function_value( &::SireVol::BoxPatching::patchIndex );
            
            BoxPatching_exposer.def( 
                "patchIndex"
                , patchIndex_function_value
                , ( bp::arg("point") )
                , "Return the index of the patch that contains the passed point" );
        
        }
        { //::SireVol::BoxPatching::patchIndexAndCenter
        
            typedef ::QPair< int, SireMaths::Vector > ( ::SireVol::BoxPatching::*patchIndexAndCenter_function_type)( ::SireMaths::Vector const & ) const;
            patchIndexAndCenter_function_type patchIndexAndCenter_function_value( &::SireVol::BoxPatching::patchIndexAndCenter );
            
            BoxPatching_exposer.def( 
                "patchIndexAndCenter"
                , patchIndexAndCenter_function_value
                , ( bp::arg("point") )
                , "Return the index of the patch that contains the passed point, together\nwith the point mapped into that patch (if the space is periodic)" );
        
        }
        { //::SireVol::BoxPatching::patchSize
        
            typedef ::SireUnits::Dimension::Length ( ::SireVol::BoxPatching::*patchSize_function_type)(  ) const;
            patchSize_function_type patchSize_function_value( &::SireVol::BoxPatching::patchSize );
            
            BoxPatching_exposer.def( 
                "patchSize"
                , patchSize_function_value
                , "Return the patch size used to control the rough size of the patches" );
        
        }
        { //::SireVol::BoxPatching::repatch
        
            typedef ::SireVol::PatchingPtr ( ::SireVol::BoxPatching::*repatch_function_type)( ::SireVol::Space const & ) const;
            repatch_function_type repatch_function_value( &::SireVol::BoxPatching::repatch );
            
            BoxPatching_exposer.def( 
                "repatch"
                , repatch_function_value
                , ( bp::arg("new_space") )
                , "Recreate the patching for the passed space" );
        
        }
        { //::SireVol::BoxPatching::toString
        
            typedef ::QString ( ::SireVol::BoxPatching::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireVol::BoxPatching::toString );
            
            BoxPatching_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of the patching" );
        
        }
        { //::SireVol::BoxPatching::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireVol::BoxPatching::typeName );
            
            BoxPatching_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        BoxPatching_exposer.staticmethod( "typeName" );
        BoxPatching_exposer.def( "__copy__", &__copy__);
        BoxPatching_exposer.def( "__deepcopy__", &__copy__);
        BoxPatching_exposer.def( "clone", &__copy__);
        BoxPatching_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireVol::BoxPatching >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BoxPatching_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireVol::BoxPatching >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BoxPatching_exposer.def( "__setstate__", &__setstate__base64< ::SireVol::BoxPatching > );
        BoxPatching_exposer.def( "__getstate__", &__getstate__base64< ::SireVol::BoxPatching > );
        BoxPatching_exposer.def( "__str__", &__str__< ::SireVol::BoxPatching > );
        BoxPatching_exposer.def( "__repr__", &__str__< ::SireVol::BoxPatching > );
    }

}
