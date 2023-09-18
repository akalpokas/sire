// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "BoreschRestraints.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "boreschrestraints.h"

#include <QDebug>

#include "boreschrestraints.h"

SireMM::BoreschRestraints __copy__(const SireMM::BoreschRestraints &other){ return SireMM::BoreschRestraints(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_BoreschRestraints_class(){

    { //::SireMM::BoreschRestraints
        typedef bp::class_< SireMM::BoreschRestraints, bp::bases< SireMM::Restraints, SireBase::Property > > BoreschRestraints_exposer_t;
        BoreschRestraints_exposer_t BoreschRestraints_exposer = BoreschRestraints_exposer_t( "BoreschRestraints", "This class provides the information for a collection of positional\nrestraints that can be added to a collection of molecues. Each\nrestraint can act on a particle or the centroid of a collection\nof particles. The restaints are spherically symmetric, and\nare either flat-bottom harmonics or harmonic potentials\n", bp::init< >("Null constructor") );
        bp::scope BoreschRestraints_scope( BoreschRestraints_exposer );
        BoreschRestraints_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        BoreschRestraints_exposer.def( bp::init< SireMM::BoreschRestraint const & >(( bp::arg("restraint") ), "") );
        BoreschRestraints_exposer.def( bp::init< QList< SireMM::BoreschRestraint > const & >(( bp::arg("restraints") ), "") );
        BoreschRestraints_exposer.def( bp::init< QString const &, SireMM::BoreschRestraint const & >(( bp::arg("name"), bp::arg("restraint") ), "") );
        BoreschRestraints_exposer.def( bp::init< QString const &, QList< SireMM::BoreschRestraint > const & >(( bp::arg("name"), bp::arg("restraints") ), "") );
        BoreschRestraints_exposer.def( bp::init< SireMM::BoreschRestraints const & >(( bp::arg("other") ), "") );
        { //::SireMM::BoreschRestraints::add
        
            typedef void ( ::SireMM::BoreschRestraints::*add_function_type)( ::SireMM::BoreschRestraint const & ) ;
            add_function_type add_function_value( &::SireMM::BoreschRestraints::add );
            
            BoreschRestraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraint") )
                , bp::release_gil_policy()
                , "Add a restraint onto the list" );
        
        }
        { //::SireMM::BoreschRestraints::add
        
            typedef void ( ::SireMM::BoreschRestraints::*add_function_type)( ::SireMM::BoreschRestraints const & ) ;
            add_function_type add_function_value( &::SireMM::BoreschRestraints::add );
            
            BoreschRestraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraints") )
                , bp::release_gil_policy()
                , "Add a restraint onto the list" );
        
        }
        { //::SireMM::BoreschRestraints::at
        
            typedef ::SireMM::BoreschRestraint const & ( ::SireMM::BoreschRestraints::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMM::BoreschRestraints::at );
            
            BoreschRestraints_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the ith restraint" );
        
        }
        { //::SireMM::BoreschRestraints::count
        
            typedef int ( ::SireMM::BoreschRestraints::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::BoreschRestraints::count );
            
            BoreschRestraints_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        { //::SireMM::BoreschRestraints::isEmpty
        
            typedef bool ( ::SireMM::BoreschRestraints::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::BoreschRestraints::isEmpty );
            
            BoreschRestraints_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMM::BoreschRestraints::isNull
        
            typedef bool ( ::SireMM::BoreschRestraints::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMM::BoreschRestraints::isNull );
            
            BoreschRestraints_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMM::BoreschRestraints::nRestraints
        
            typedef int ( ::SireMM::BoreschRestraints::*nRestraints_function_type)(  ) const;
            nRestraints_function_type nRestraints_function_value( &::SireMM::BoreschRestraints::nRestraints );
            
            BoreschRestraints_exposer.def( 
                "nRestraints"
                , nRestraints_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        BoreschRestraints_exposer.def( bp::self != bp::self );
        BoreschRestraints_exposer.def( bp::self + bp::other< SireMM::BoreschRestraint >() );
        BoreschRestraints_exposer.def( bp::self + bp::self );
        { //::SireMM::BoreschRestraints::operator=
        
            typedef ::SireMM::BoreschRestraints & ( ::SireMM::BoreschRestraints::*assign_function_type)( ::SireMM::BoreschRestraints const & ) ;
            assign_function_type assign_function_value( &::SireMM::BoreschRestraints::operator= );
            
            BoreschRestraints_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        BoreschRestraints_exposer.def( bp::self == bp::self );
        { //::SireMM::BoreschRestraints::operator[]
        
            typedef ::SireMM::BoreschRestraint const & ( ::SireMM::BoreschRestraints::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::BoreschRestraints::operator[] );
            
            BoreschRestraints_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::BoreschRestraints::restraints
        
            typedef ::QList< SireMM::BoreschRestraint > ( ::SireMM::BoreschRestraints::*restraints_function_type)(  ) const;
            restraints_function_type restraints_function_value( &::SireMM::BoreschRestraints::restraints );
            
            BoreschRestraints_exposer.def( 
                "restraints"
                , restraints_function_value
                , bp::release_gil_policy()
                , "Return all of the restraints" );
        
        }
        { //::SireMM::BoreschRestraints::size
        
            typedef int ( ::SireMM::BoreschRestraints::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::BoreschRestraints::size );
            
            BoreschRestraints_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        { //::SireMM::BoreschRestraints::toString
        
            typedef ::QString ( ::SireMM::BoreschRestraints::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::BoreschRestraints::toString );
            
            BoreschRestraints_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::BoreschRestraints::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::BoreschRestraints::typeName );
            
            BoreschRestraints_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::BoreschRestraints::what
        
            typedef char const * ( ::SireMM::BoreschRestraints::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::BoreschRestraints::what );
            
            BoreschRestraints_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        BoreschRestraints_exposer.staticmethod( "typeName" );
        BoreschRestraints_exposer.def( "__copy__", &__copy__);
        BoreschRestraints_exposer.def( "__deepcopy__", &__copy__);
        BoreschRestraints_exposer.def( "clone", &__copy__);
        BoreschRestraints_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::BoreschRestraints >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BoreschRestraints_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::BoreschRestraints >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BoreschRestraints_exposer.def_pickle(sire_pickle_suite< ::SireMM::BoreschRestraints >());
        BoreschRestraints_exposer.def( "__str__", &__str__< ::SireMM::BoreschRestraints > );
        BoreschRestraints_exposer.def( "__repr__", &__str__< ::SireMM::BoreschRestraints > );
        BoreschRestraints_exposer.def( "__len__", &__len_size< ::SireMM::BoreschRestraints > );
    }

}
