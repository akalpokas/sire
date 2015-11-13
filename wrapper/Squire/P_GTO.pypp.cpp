// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "P_GTO.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include "SireMaths/boys.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "pgto.h"

#include "pointcharge.h"

#include "pointdipole.h"

#include "pgto.h"

Squire::P_GTO __copy__(const Squire::P_GTO &other){ return Squire::P_GTO(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_P_GTO_class(){

    { //::Squire::P_GTO
        typedef bp::class_< Squire::P_GTO, bp::bases< Squire::GTO, Squire::OrbitalShell, Squire::Orbital, SireBase::Property > > P_GTO_exposer_t;
        P_GTO_exposer_t P_GTO_exposer = P_GTO_exposer_t( "P_GTO", bp::init< >() );
        bp::scope P_GTO_scope( P_GTO_exposer );
        P_GTO_exposer.def( bp::init< double, bp::optional< double > >(( bp::arg("alpha"), bp::arg("scale")=1 )) );
        P_GTO_exposer.def( bp::init< Squire::P_GTO const & >(( bp::arg("other") )) );
        { //::Squire::P_GTO::angularMomentum
        
            typedef int ( ::Squire::P_GTO::*angularMomentum_function_type)(  ) const;
            angularMomentum_function_type angularMomentum_function_value( &::Squire::P_GTO::angularMomentum );
            
            P_GTO_exposer.def( 
                "angularMomentum"
                , angularMomentum_function_value );
        
        }
        { //::Squire::P_GTO::nOrbitals
        
            typedef int ( ::Squire::P_GTO::*nOrbitals_function_type)(  ) const;
            nOrbitals_function_type nOrbitals_function_value( &::Squire::P_GTO::nOrbitals );
            
            P_GTO_exposer.def( 
                "nOrbitals"
                , nOrbitals_function_value );
        
        }
        P_GTO_exposer.def( bp::self != bp::self );
        { //::Squire::P_GTO::operator=
        
            typedef ::Squire::P_GTO & ( ::Squire::P_GTO::*assign_function_type)( ::Squire::P_GTO const & ) ;
            assign_function_type assign_function_value( &::Squire::P_GTO::operator= );
            
            P_GTO_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        P_GTO_exposer.def( bp::self == bp::self );
        { //::Squire::P_GTO::toString
        
            typedef ::QString ( ::Squire::P_GTO::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::Squire::P_GTO::toString );
            
            P_GTO_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::Squire::P_GTO::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::P_GTO::typeName );
            
            P_GTO_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        P_GTO_exposer.staticmethod( "typeName" );
        P_GTO_exposer.def( "__copy__", &__copy__);
        P_GTO_exposer.def( "__deepcopy__", &__copy__);
        P_GTO_exposer.def( "clone", &__copy__);
        P_GTO_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::P_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        P_GTO_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::P_GTO >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        P_GTO_exposer.def( "__str__", &__str__< ::Squire::P_GTO > );
        P_GTO_exposer.def( "__repr__", &__str__< ::Squire::P_GTO > );
    }

}
