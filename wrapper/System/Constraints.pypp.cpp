// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Constraints.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/savestate.h"

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMol/molecules.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/errors.h"

#include "constraints.h"

#include "delta.h"

#include "moleculeconstraint.h"

#include "system.h"

#include <QDebug>

#include "constraints.h"

SireSystem::Constraints __copy__(const SireSystem::Constraints &other){ return SireSystem::Constraints(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Constraints_class(){

    { //::SireSystem::Constraints
        typedef bp::class_< SireSystem::Constraints, bp::bases< SireBase::Property > > Constraints_exposer_t;
        Constraints_exposer_t Constraints_exposer = Constraints_exposer_t( "Constraints", "This class is used to hold and organise a collection of Constraint objects\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope Constraints_scope( Constraints_exposer );
        Constraints_exposer.def( bp::init< SireSystem::Constraint const & >(( bp::arg("constraint") ), "Construct to contain just the constraint constraint") );
        Constraints_exposer.def( bp::init< QVector< SireBase::PropPtr< SireSystem::Constraint > > const & >(( bp::arg("constraints") ), "Construct from the passed list of constraints") );
        Constraints_exposer.def( bp::init< QList< SireBase::PropPtr< SireSystem::Constraint > > const & >(( bp::arg("constraints") ), "Construct from the passed list of constraints") );
        Constraints_exposer.def( bp::init< SireSystem::Constraints const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::Constraints::add
        
            typedef void ( ::SireSystem::Constraints::*add_function_type)( ::SireSystem::Constraint const & ) ;
            add_function_type add_function_value( &::SireSystem::Constraints::add );
            
            Constraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("constraint") )
                , bp::release_gil_policy()
                , "Add the passed constraint to this set - this is only added\nif it does not exist in this set" );
        
        }
        { //::SireSystem::Constraints::add
        
            typedef void ( ::SireSystem::Constraints::*add_function_type)( ::SireSystem::Constraints const & ) ;
            add_function_type add_function_value( &::SireSystem::Constraints::add );
            
            Constraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("constraints") )
                , bp::release_gil_policy()
                , "Add all of the passed constraints to this set. This only\nadds the constraints that are not already part of this set" );
        
        }
        { //::SireSystem::Constraints::apply
        
            typedef ::SireSystem::System ( ::SireSystem::Constraints::*apply_function_type)( ::SireSystem::System const & ) ;
            apply_function_type apply_function_value( &::SireSystem::Constraints::apply );
            
            Constraints_exposer.def( 
                "apply"
                , apply_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Apply all of the constraints to the passed system. This\nreturns a system that satisfies all of the constraints" );
        
        }
        { //::SireSystem::Constraints::areSatisfied
        
            typedef bool ( ::SireSystem::Constraints::*areSatisfied_function_type)( ::SireSystem::System const & ) const;
            areSatisfied_function_type areSatisfied_function_value( &::SireSystem::Constraints::areSatisfied );
            
            Constraints_exposer.def( 
                "areSatisfied"
                , areSatisfied_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Return whether or not all of the constraints in this set are\nsatisfied in the passed system" );
        
        }
        { //::SireSystem::Constraints::assertSatisfied
        
            typedef void ( ::SireSystem::Constraints::*assertSatisfied_function_type)( ::SireSystem::System const & ) const;
            assertSatisfied_function_type assertSatisfied_function_value( &::SireSystem::Constraints::assertSatisfied );
            
            Constraints_exposer.def( 
                "assertSatisfied"
                , assertSatisfied_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Assert that all of the constraints in this set are satisfied\nin the passed system\nThrow: SireSystem::constraint_error\n" );
        
        }
        { //::SireSystem::Constraints::constraints
        
            typedef ::QVector< SireBase::PropPtr< SireSystem::Constraint > > ( ::SireSystem::Constraints::*constraints_function_type)(  ) const;
            constraints_function_type constraints_function_value( &::SireSystem::Constraints::constraints );
            
            Constraints_exposer.def( 
                "constraints"
                , constraints_function_value
                , bp::release_gil_policy()
                , "Return the list of all of the constraints in this set" );
        
        }
        { //::SireSystem::Constraints::count
        
            typedef int ( ::SireSystem::Constraints::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireSystem::Constraints::count );
            
            Constraints_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of constraints in this set" );
        
        }
        { //::SireSystem::Constraints::isEmpty
        
            typedef bool ( ::SireSystem::Constraints::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireSystem::Constraints::isEmpty );
            
            Constraints_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether this is empty (contains no constraints)" );
        
        }
        { //::SireSystem::Constraints::nConstraints
        
            typedef int ( ::SireSystem::Constraints::*nConstraints_function_type)(  ) const;
            nConstraints_function_type nConstraints_function_value( &::SireSystem::Constraints::nConstraints );
            
            Constraints_exposer.def( 
                "nConstraints"
                , nConstraints_function_value
                , bp::release_gil_policy()
                , "Return the number of constraints in this set" );
        
        }
        Constraints_exposer.def( bp::self != bp::self );
        { //::SireSystem::Constraints::operator=
        
            typedef ::SireSystem::Constraints & ( ::SireSystem::Constraints::*assign_function_type)( ::SireSystem::Constraints const & ) ;
            assign_function_type assign_function_value( &::SireSystem::Constraints::operator= );
            
            Constraints_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Constraints_exposer.def( bp::self == bp::self );
        { //::SireSystem::Constraints::operator[]
        
            typedef ::SireSystem::Constraint const & ( ::SireSystem::Constraints::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireSystem::Constraints::operator[] );
            
            Constraints_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireSystem::Constraints::remove
        
            typedef void ( ::SireSystem::Constraints::*remove_function_type)( ::SireSystem::Constraint const & ) ;
            remove_function_type remove_function_value( &::SireSystem::Constraints::remove );
            
            Constraints_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("constraint") )
                , bp::release_gil_policy()
                , "Remove the constraint constraint from this set - this\ndoes nothing if this constraint is not part of this set" );
        
        }
        { //::SireSystem::Constraints::remove
        
            typedef void ( ::SireSystem::Constraints::*remove_function_type)( ::SireSystem::Constraints const & ) ;
            remove_function_type remove_function_value( &::SireSystem::Constraints::remove );
            
            Constraints_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("constraints") )
                , bp::release_gil_policy()
                , "Remove all of the constraints in constraints from this\nset - this ignores constraints that are not in this set" );
        
        }
        { //::SireSystem::Constraints::removeAt
        
            typedef void ( ::SireSystem::Constraints::*removeAt_function_type)( int ) ;
            removeAt_function_type removeAt_function_value( &::SireSystem::Constraints::removeAt );
            
            Constraints_exposer.def( 
                "removeAt"
                , removeAt_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Remove the constraint at index i" );
        
        }
        { //::SireSystem::Constraints::size
        
            typedef int ( ::SireSystem::Constraints::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireSystem::Constraints::size );
            
            Constraints_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of constraints in this set" );
        
        }
        { //::SireSystem::Constraints::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::Constraints::typeName );
            
            Constraints_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireSystem::Constraints::what
        
            typedef char const * ( ::SireSystem::Constraints::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireSystem::Constraints::what );
            
            Constraints_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Constraints_exposer.staticmethod( "typeName" );
        Constraints_exposer.def( "__copy__", &__copy__<SireSystem::Constraints>);
        Constraints_exposer.def( "__deepcopy__", &__copy__<SireSystem::Constraints>);
        Constraints_exposer.def( "clone", &__copy__<SireSystem::Constraints>);
        Constraints_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::Constraints >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Constraints_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::Constraints >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Constraints_exposer.def_pickle(sire_pickle_suite< ::SireSystem::Constraints >());
        Constraints_exposer.def( "__str__", &__str__< ::SireSystem::Constraints > );
        Constraints_exposer.def( "__repr__", &__str__< ::SireSystem::Constraints > );
        Constraints_exposer.def( "__len__", &__len_size< ::SireSystem::Constraints > );
    }

}
