// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CheckPoint.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "checkpoint.h"

#include "checkpoint.h"

SireSystem::CheckPoint __copy__(const SireSystem::CheckPoint &other){ return SireSystem::CheckPoint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CheckPoint_class(){

    { //::SireSystem::CheckPoint
        typedef bp::class_< SireSystem::CheckPoint, bp::bases< SireBase::Property > > CheckPoint_exposer_t;
        CheckPoint_exposer_t CheckPoint_exposer = CheckPoint_exposer_t( "CheckPoint", "This class holds a checkpoint of a system. This allows you to\nsave the current state of a system so that you can restore\nit at a later point (or even save it to diskdatabase)\n\nAuthor: Christopher Woods\n", bp::init< >("Create a null (empty) checkpoint") );
        bp::scope CheckPoint_scope( CheckPoint_exposer );
        CheckPoint_exposer.def( bp::init< SireSystem::System const & >(( bp::arg("system") ), "Construct a checkpoint to hold the current state of the system system") );
        CheckPoint_exposer.def( bp::init< SireSystem::CheckPoint const & >(( bp::arg("other") ), "Copy constructor") );
        CheckPoint_exposer.def( bp::self != bp::self );
        { //::SireSystem::CheckPoint::operator=
        
            typedef ::SireSystem::CheckPoint & ( ::SireSystem::CheckPoint::*assign_function_type)( ::SireSystem::System const & ) ;
            assign_function_type assign_function_value( &::SireSystem::CheckPoint::operator= );
            
            CheckPoint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("system") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireSystem::CheckPoint::operator=
        
            typedef ::SireSystem::CheckPoint & ( ::SireSystem::CheckPoint::*assign_function_type)( ::SireSystem::CheckPoint const & ) ;
            assign_function_type assign_function_value( &::SireSystem::CheckPoint::operator= );
            
            CheckPoint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CheckPoint_exposer.def( bp::self == bp::self );
        { //::SireSystem::CheckPoint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::CheckPoint::typeName );
            
            CheckPoint_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        CheckPoint_exposer.staticmethod( "typeName" );
        CheckPoint_exposer.def( "__copy__", &__copy__);
        CheckPoint_exposer.def( "__deepcopy__", &__copy__);
        CheckPoint_exposer.def( "clone", &__copy__);
        CheckPoint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::CheckPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CheckPoint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::CheckPoint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CheckPoint_exposer.def( "__setstate__", &__setstate__base64< ::SireSystem::CheckPoint > );
        CheckPoint_exposer.def( "__getstate__", &__getstate__base64< ::SireSystem::CheckPoint > );
        CheckPoint_exposer.def( "__str__", &__str__< ::SireSystem::CheckPoint > );
        CheckPoint_exposer.def( "__repr__", &__str__< ::SireSystem::CheckPoint > );
    }

}
