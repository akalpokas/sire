// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "SysID.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireSystem/errors.h"

#include "sysid.h"

#include "sysidx.h"

#include "sysname.h"

#include "sysid.h"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_SysID_class(){

    { //::SireSystem::SysID
        typedef bp::class_< SireSystem::SysID, bp::bases< SireID::ID >, boost::noncopyable > SysID_exposer_t;
        SysID_exposer_t SysID_exposer = SysID_exposer_t( "SysID", "The base class of all system identifiers\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope SysID_scope( SysID_exposer );
        { //::SireSystem::SysID::map
        
            typedef ::QList< SireSystem::SysIdx > ( ::SireSystem::SysID::*map_function_type)( ::SireSystem::Systems const & ) const;
            map_function_type map_function_value( &::SireSystem::SysID::map );
            
            SysID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("systems") )
                , bp::release_gil_policy()
                , "" );
        
        }
        SysID_exposer.def( bp::self & bp::self );
        { //::SireSystem::SysID::operator()
        
            typedef ::SireID::Specify< SireSystem::SysID > ( ::SireSystem::SysID::*__call___function_type)( int ) const;
            __call___function_type __call___function_value( &::SireSystem::SysID::operator() );
            
            SysID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireSystem::SysID::operator()
        
            typedef ::SireID::Specify< SireSystem::SysID > ( ::SireSystem::SysID::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireSystem::SysID::operator() );
            
            SysID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , "" );
        
        }
        SysID_exposer.def( bp::self * bp::self );
        SysID_exposer.def( bp::self + bp::self );
        { //::SireSystem::SysID::operator[]
        
            typedef ::SireID::Specify< SireSystem::SysID > ( ::SireSystem::SysID::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireSystem::SysID::operator[] );
            
            SysID_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        SysID_exposer.def( bp::self | bp::self );
        { //::SireSystem::SysID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::SysID::typeName );
            
            SysID_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SysID_exposer.staticmethod( "typeName" );
        SysID_exposer.def( "__str__", &__str__< ::SireSystem::SysID > );
        SysID_exposer.def( "__repr__", &__str__< ::SireSystem::SysID > );
        SysID_exposer.def( "__hash__", &::SireSystem::SysID::hash );
    }

}
