// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "AngleRestraints.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/conditional.h"

#include "SireCAS/errors.h"

#include "SireCAS/power.h"

#include "SireCAS/symbols.h"

#include "SireCAS/values.h"

#include "SireFF/forcetable.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/angle.h"

#include "SireUnits/units.h"

#include "anglerestraint.h"

#include <QDebug>

#include "anglerestraint.h"

SireMM::AngleRestraints __copy__(const SireMM::AngleRestraints &other){ return SireMM::AngleRestraints(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_AngleRestraints_class(){

    { //::SireMM::AngleRestraints
        typedef bp::class_< SireMM::AngleRestraints, bp::bases< SireMM::Restraints, SireBase::Property > > AngleRestraints_exposer_t;
        AngleRestraints_exposer_t AngleRestraints_exposer = AngleRestraints_exposer_t( "AngleRestraints", "This class represents a collection of angle restraints", bp::init< >("Null constructor") );
        bp::scope AngleRestraints_scope( AngleRestraints_exposer );
        AngleRestraints_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        AngleRestraints_exposer.def( bp::init< SireMM::AngleRestraint const & >(( bp::arg("restraint") ), "") );
        AngleRestraints_exposer.def( bp::init< QList< SireMM::AngleRestraint > const & >(( bp::arg("restraints") ), "") );
        AngleRestraints_exposer.def( bp::init< QString const &, SireMM::AngleRestraint const & >(( bp::arg("name"), bp::arg("restraint") ), "") );
        AngleRestraints_exposer.def( bp::init< QString const &, QList< SireMM::AngleRestraint > const & >(( bp::arg("name"), bp::arg("restraints") ), "") );
        AngleRestraints_exposer.def( bp::init< SireMM::AngleRestraints const & >(( bp::arg("other") ), "") );
        { //::SireMM::AngleRestraints::add
        
            typedef void ( ::SireMM::AngleRestraints::*add_function_type)( ::SireMM::AngleRestraint const & ) ;
            add_function_type add_function_value( &::SireMM::AngleRestraints::add );
            
            AngleRestraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraint") )
                , bp::release_gil_policy()
                , "Add a restraint onto the list" );
        
        }
        { //::SireMM::AngleRestraints::add
        
            typedef void ( ::SireMM::AngleRestraints::*add_function_type)( ::SireMM::AngleRestraints const & ) ;
            add_function_type add_function_value( &::SireMM::AngleRestraints::add );
            
            AngleRestraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraints") )
                , bp::release_gil_policy()
                , "Add a restraint onto the list" );
        
        }
        { //::SireMM::AngleRestraints::at
        
            typedef ::SireMM::AngleRestraint const & ( ::SireMM::AngleRestraints::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMM::AngleRestraints::at );
            
            AngleRestraints_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the ith restraint" );
        
        }
        { //::SireMM::AngleRestraints::count
        
            typedef int ( ::SireMM::AngleRestraints::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::AngleRestraints::count );
            
            AngleRestraints_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        { //::SireMM::AngleRestraints::isEmpty
        
            typedef bool ( ::SireMM::AngleRestraints::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::AngleRestraints::isEmpty );
            
            AngleRestraints_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMM::AngleRestraints::isNull
        
            typedef bool ( ::SireMM::AngleRestraints::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMM::AngleRestraints::isNull );
            
            AngleRestraints_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMM::AngleRestraints::nRestraints
        
            typedef int ( ::SireMM::AngleRestraints::*nRestraints_function_type)(  ) const;
            nRestraints_function_type nRestraints_function_value( &::SireMM::AngleRestraints::nRestraints );
            
            AngleRestraints_exposer.def( 
                "nRestraints"
                , nRestraints_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        AngleRestraints_exposer.def( bp::self != bp::self );
        AngleRestraints_exposer.def( bp::self + bp::other< SireMM::AngleRestraint >() );
        AngleRestraints_exposer.def( bp::self + bp::self );
        { //::SireMM::AngleRestraints::operator=
        
            typedef ::SireMM::AngleRestraints & ( ::SireMM::AngleRestraints::*assign_function_type)( ::SireMM::AngleRestraints const & ) ;
            assign_function_type assign_function_value( &::SireMM::AngleRestraints::operator= );
            
            AngleRestraints_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AngleRestraints_exposer.def( bp::self == bp::self );
        { //::SireMM::AngleRestraints::operator[]
        
            typedef ::SireMM::AngleRestraint const & ( ::SireMM::AngleRestraints::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::AngleRestraints::operator[] );
            
            AngleRestraints_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::AngleRestraints::restraints
        
            typedef ::QList< SireMM::AngleRestraint > ( ::SireMM::AngleRestraints::*restraints_function_type)(  ) const;
            restraints_function_type restraints_function_value( &::SireMM::AngleRestraints::restraints );
            
            AngleRestraints_exposer.def( 
                "restraints"
                , restraints_function_value
                , bp::release_gil_policy()
                , "Return all of the restraints" );
        
        }
        { //::SireMM::AngleRestraints::size
        
            typedef int ( ::SireMM::AngleRestraints::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::AngleRestraints::size );
            
            AngleRestraints_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        { //::SireMM::AngleRestraints::toString
        
            typedef ::QString ( ::SireMM::AngleRestraints::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::AngleRestraints::toString );
            
            AngleRestraints_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::AngleRestraints::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::AngleRestraints::typeName );
            
            AngleRestraints_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::AngleRestraints::what
        
            typedef char const * ( ::SireMM::AngleRestraints::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::AngleRestraints::what );
            
            AngleRestraints_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AngleRestraints_exposer.staticmethod( "typeName" );
        AngleRestraints_exposer.def( "__copy__", &__copy__<SireMM::AngleRestraints>);
        AngleRestraints_exposer.def( "__deepcopy__", &__copy__<SireMM::AngleRestraints>);
        AngleRestraints_exposer.def( "clone", &__copy__<SireMM::AngleRestraints>);
        AngleRestraints_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::AngleRestraints >,
                                bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AngleRestraints_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::AngleRestraints >,
                                bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AngleRestraints_exposer.def_pickle(sire_pickle_suite< ::SireMM::AngleRestraints >());
        AngleRestraints_exposer.def( "__str__", &__str__< ::SireMM::AngleRestraints > );
        AngleRestraints_exposer.def( "__repr__", &__str__< ::SireMM::AngleRestraints > );
        AngleRestraints_exposer.def( "__len__", &__len_size< ::SireMM::AngleRestraints > );
    }

}
