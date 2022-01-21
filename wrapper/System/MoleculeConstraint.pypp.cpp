// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MoleculeConstraint.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/molecules.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "delta.h"

#include "moleculeconstraint.h"

#include "system.h"

#include <QDebug>

#include <boost/shared_ptr.hpp>

#include "moleculeconstraint.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_MoleculeConstraint_class(){

    { //::SireSystem::MoleculeConstraint
        typedef bp::class_< SireSystem::MoleculeConstraint, bp::bases< SireSystem::Constraint, SireBase::Property >, boost::noncopyable > MoleculeConstraint_exposer_t;
        MoleculeConstraint_exposer_t MoleculeConstraint_exposer = MoleculeConstraint_exposer_t( "MoleculeConstraint", "This is the base class of all Molecule constraints. These\nare constraints that affect molecules in a system. A molecule\nconstraint works by being updated by a system and returning\nthe molecules that must then be changed to maintain the\nconstraint\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope MoleculeConstraint_scope( MoleculeConstraint_exposer );
        { //::SireSystem::MoleculeConstraint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::MoleculeConstraint::typeName );
            
            MoleculeConstraint_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        MoleculeConstraint_exposer.staticmethod( "typeName" );
        MoleculeConstraint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::MoleculeConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MoleculeConstraint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::MoleculeConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MoleculeConstraint_exposer.def( "__setstate__", &__setstate__base64< ::SireSystem::MoleculeConstraint > );
        MoleculeConstraint_exposer.def( "__getstate__", &__getstate__base64< ::SireSystem::MoleculeConstraint > );
        MoleculeConstraint_exposer.def( "__str__", &__str__< ::SireSystem::MoleculeConstraint > );
        MoleculeConstraint_exposer.def( "__repr__", &__str__< ::SireSystem::MoleculeConstraint > );
    }

}
