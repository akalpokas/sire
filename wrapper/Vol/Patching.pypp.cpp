// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Patching.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "patching.h"

#include "patching.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Patching_class(){

    { //::SireVol::Patching
        typedef bp::class_< SireVol::Patching, bp::bases< SireBase::Property >, boost::noncopyable > Patching_exposer_t;
        Patching_exposer_t Patching_exposer = Patching_exposer_t( "Patching", "This is the base class of all Patching classes. Patching\nrepresents a scheme for decomposing a space into a set\nof regions (domain decomposition) that contain neighbouring\nCoordGroups - this allows inter-CoordGroup calculations to\nbe accelerated as cutoff tests can be applied to the patches\nto eliminate tests of the contained CoordGroups.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Patching_scope( Patching_exposer );
        { //::SireVol::Patching::nPatches
        
            typedef int ( ::SireVol::Patching::*nPatches_function_type)(  ) const;
            nPatches_function_type nPatches_function_value( &::SireVol::Patching::nPatches );
            
            Patching_exposer.def( 
                "nPatches"
                , nPatches_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Patching::null
        
            typedef ::SireVol::NullPatching ( *null_function_type )(  );
            null_function_type null_function_value( &::SireVol::Patching::null );
            
            Patching_exposer.def( 
                "null"
                , null_function_value
                , bp::release_gil_policy()
                , "Retunr the null patching object" );
        
        }
        { //::SireVol::Patching::patchIndex
        
            typedef int ( ::SireVol::Patching::*patchIndex_function_type)( ::SireMaths::Vector const & ) const;
            patchIndex_function_type patchIndex_function_value( &::SireVol::Patching::patchIndex );
            
            Patching_exposer.def( 
                "patchIndex"
                , patchIndex_function_value
                , ( bp::arg("point") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Patching::patchIndexAndCenter
        
            typedef ::QPair< int, SireMaths::Vector > ( ::SireVol::Patching::*patchIndexAndCenter_function_type)( ::SireMaths::Vector const & ) const;
            patchIndexAndCenter_function_type patchIndexAndCenter_function_value( &::SireVol::Patching::patchIndexAndCenter );
            
            Patching_exposer.def( 
                "patchIndexAndCenter"
                , patchIndexAndCenter_function_value
                , ( bp::arg("point") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Patching::rebalance
        
            typedef ::SireVol::PatchingPtr ( ::SireVol::Patching::*rebalance_function_type)( ::SireVol::Space const &,::QVector< SireVol::CoordGroupArray > const & ) const;
            rebalance_function_type rebalance_function_value( &::SireVol::Patching::rebalance );
            
            Patching_exposer.def( 
                "rebalance"
                , rebalance_function_value
                , ( bp::arg("space"), bp::arg("patchcoords") )
                , bp::release_gil_policy()
                , "Rebalance the patching so that the patches for the passed space contain\nroughly equal numbers of CoordGroups" );
        
        }
        { //::SireVol::Patching::repatch
        
            typedef ::SireVol::PatchingPtr ( ::SireVol::Patching::*repatch_function_type)( ::SireVol::Space const & ) const;
            repatch_function_type repatch_function_value( &::SireVol::Patching::repatch );
            
            Patching_exposer.def( 
                "repatch"
                , repatch_function_value
                , ( bp::arg("new_space") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::Patching::space
        
            typedef ::SireVol::Space const & ( ::SireVol::Patching::*space_function_type)(  ) const;
            space_function_type space_function_value( &::SireVol::Patching::space );
            
            Patching_exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the space used to create this patching scheme" );
        
        }
        { //::SireVol::Patching::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireVol::Patching::typeName );
            
            Patching_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Patching_exposer.staticmethod( "null" );
        Patching_exposer.staticmethod( "typeName" );
        Patching_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireVol::Patching >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Patching_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireVol::Patching >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Patching_exposer.def_pickle(sire_pickle_suite< ::SireVol::Patching >());
        Patching_exposer.def( "__str__", &__str__< ::SireVol::Patching > );
        Patching_exposer.def( "__repr__", &__str__< ::SireVol::Patching > );
    }

}
