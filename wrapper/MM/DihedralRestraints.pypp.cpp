// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "DihedralRestraints.pypp.hpp"

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

#include "dihedralrestraint.h"

#include <QDebug>

#include "dihedralrestraint.h"

SireMM::DihedralRestraints __copy__(const SireMM::DihedralRestraints &other){ return SireMM::DihedralRestraints(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_DihedralRestraints_class(){

    { //::SireMM::DihedralRestraints
        typedef bp::class_< SireMM::DihedralRestraints, bp::bases< SireMM::Restraints, SireBase::Property > > DihedralRestraints_exposer_t;
        DihedralRestraints_exposer_t DihedralRestraints_exposer = DihedralRestraints_exposer_t( "DihedralRestraints", "This class represents a collection of angle restraints", bp::init< >("Null constructor") );
        bp::scope DihedralRestraints_scope( DihedralRestraints_exposer );
        DihedralRestraints_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        DihedralRestraints_exposer.def( bp::init< SireMM::DihedralRestraint const & >(( bp::arg("restraint") ), "") );
        DihedralRestraints_exposer.def( bp::init< QList< SireMM::DihedralRestraint > const & >(( bp::arg("restraints") ), "") );
        DihedralRestraints_exposer.def( bp::init< QString const &, SireMM::DihedralRestraint const & >(( bp::arg("name"), bp::arg("restraint") ), "") );
        DihedralRestraints_exposer.def( bp::init< QString const &, QList< SireMM::DihedralRestraint > const & >(( bp::arg("name"), bp::arg("restraints") ), "") );
        DihedralRestraints_exposer.def( bp::init< SireMM::DihedralRestraints const & >(( bp::arg("other") ), "") );
        { //::SireMM::DihedralRestraints::add
        
            typedef void ( ::SireMM::DihedralRestraints::*add_function_type)( ::SireMM::DihedralRestraint const & ) ;
            add_function_type add_function_value( &::SireMM::DihedralRestraints::add );
            
            DihedralRestraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraint") )
                , bp::release_gil_policy()
                , "Add a restraints onto the list" );
        
        }
        { //::SireMM::DihedralRestraints::add
        
            typedef void ( ::SireMM::DihedralRestraints::*add_function_type)( ::SireMM::DihedralRestraints const & ) ;
            add_function_type add_function_value( &::SireMM::DihedralRestraints::add );
            
            DihedralRestraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraints") )
                , bp::release_gil_policy()
                , "Add a restraint onto the list" );
        
        }
        { //::SireMM::DihedralRestraints::at
        
            typedef ::SireMM::DihedralRestraint const & ( ::SireMM::DihedralRestraints::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMM::DihedralRestraints::at );
            
            DihedralRestraints_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the ith restraint" );
        
        }
        { //::SireMM::DihedralRestraints::count
        
            typedef int ( ::SireMM::DihedralRestraints::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::DihedralRestraints::count );
            
            DihedralRestraints_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        { //::SireMM::DihedralRestraints::isEmpty
        
            typedef bool ( ::SireMM::DihedralRestraints::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::DihedralRestraints::isEmpty );
            
            DihedralRestraints_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMM::DihedralRestraints::isNull
        
            typedef bool ( ::SireMM::DihedralRestraints::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMM::DihedralRestraints::isNull );
            
            DihedralRestraints_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMM::DihedralRestraints::nRestraints
        
            typedef int ( ::SireMM::DihedralRestraints::*nRestraints_function_type)(  ) const;
            nRestraints_function_type nRestraints_function_value( &::SireMM::DihedralRestraints::nRestraints );
            
            DihedralRestraints_exposer.def( 
                "nRestraints"
                , nRestraints_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        DihedralRestraints_exposer.def( bp::self != bp::self );
        DihedralRestraints_exposer.def( bp::self + bp::other< SireMM::DihedralRestraint >() );
        DihedralRestraints_exposer.def( bp::self + bp::self );
        { //::SireMM::DihedralRestraints::operator=
        
            typedef ::SireMM::DihedralRestraints & ( ::SireMM::DihedralRestraints::*assign_function_type)( ::SireMM::DihedralRestraints const & ) ;
            assign_function_type assign_function_value( &::SireMM::DihedralRestraints::operator= );
            
            DihedralRestraints_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        DihedralRestraints_exposer.def( bp::self == bp::self );
        { //::SireMM::DihedralRestraints::operator[]
        
            typedef ::SireMM::DihedralRestraint const & ( ::SireMM::DihedralRestraints::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::DihedralRestraints::operator[] );
            
            DihedralRestraints_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::DihedralRestraints::restraints
        
            typedef ::QList< SireMM::DihedralRestraint > ( ::SireMM::DihedralRestraints::*restraints_function_type)(  ) const;
            restraints_function_type restraints_function_value( &::SireMM::DihedralRestraints::restraints );
            
            DihedralRestraints_exposer.def( 
                "restraints"
                , restraints_function_value
                , bp::release_gil_policy()
                , "Return all of the restraints" );
        
        }
        { //::SireMM::DihedralRestraints::size
        
            typedef int ( ::SireMM::DihedralRestraints::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::DihedralRestraints::size );
            
            DihedralRestraints_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        { //::SireMM::DihedralRestraints::toString
        
            typedef ::QString ( ::SireMM::DihedralRestraints::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::DihedralRestraints::toString );
            
            DihedralRestraints_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::DihedralRestraints::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::DihedralRestraints::typeName );
            
            DihedralRestraints_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::DihedralRestraints::what
        
            typedef char const * ( ::SireMM::DihedralRestraints::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::DihedralRestraints::what );
            
            DihedralRestraints_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        DihedralRestraints_exposer.staticmethod( "typeName" );
        DihedralRestraints_exposer.def( "__copy__", &__copy__<SireMM::DihedralRestraints>);
        DihedralRestraints_exposer.def( "__deepcopy__", &__copy__<SireMM::DihedralRestraints>);
        DihedralRestraints_exposer.def( "clone", &__copy__<SireMM::DihedralRestraints>);
        DihedralRestraints_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::DihedralRestraints >,
                                bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DihedralRestraints_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::DihedralRestraints >,
                                bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DihedralRestraints_exposer.def_pickle(sire_pickle_suite< ::SireMM::DihedralRestraints >());
        DihedralRestraints_exposer.def( "__str__", &__str__< ::SireMM::DihedralRestraints > );
        DihedralRestraints_exposer.def( "__repr__", &__str__< ::SireMM::DihedralRestraints > );
        DihedralRestraints_exposer.def( "__len__", &__len_size< ::SireMM::DihedralRestraints > );
    }

}
