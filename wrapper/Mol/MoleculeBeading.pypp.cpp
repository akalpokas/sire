// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MoleculeBeading.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atombeads.h"

#include "atomidx.h"

#include "atomselection.h"

#include "beadidx.h"

#include "beading.h"

#include "moleculedata.h"

#include "moleculeinfodata.h"

#include <boost/noncopyable.hpp>

#include "beading.h"

SireMol::MoleculeBeading __copy__(const SireMol::MoleculeBeading &other){ return SireMol::MoleculeBeading(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_MoleculeBeading_class(){

    { //::SireMol::MoleculeBeading
        typedef bp::class_< SireMol::MoleculeBeading, bp::bases< SireMol::Beading, SireMol::MolViewProperty, SireBase::Property > > MoleculeBeading_exposer_t;
        MoleculeBeading_exposer_t MoleculeBeading_exposer = MoleculeBeading_exposer_t( "MoleculeBeading", bp::init< >() );
        bp::scope MoleculeBeading_scope( MoleculeBeading_exposer );
        MoleculeBeading_exposer.def( bp::init< SireMol::MoleculeBeading const & >(( bp::arg("other") )) );
        MoleculeBeading_exposer.def( bp::self != bp::self );
        { //::SireMol::MoleculeBeading::operator=
        
            typedef ::SireMol::MoleculeBeading & ( ::SireMol::MoleculeBeading::*assign_function_type)( ::SireMol::MoleculeBeading const & ) ;
            assign_function_type assign_function_value( &::SireMol::MoleculeBeading::operator= );
            
            MoleculeBeading_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        MoleculeBeading_exposer.def( bp::self == bp::self );
        { //::SireMol::MoleculeBeading::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::MoleculeBeading::typeName );
            
            MoleculeBeading_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        MoleculeBeading_exposer.staticmethod( "typeName" );
        MoleculeBeading_exposer.def( "__copy__", &__copy__);
        MoleculeBeading_exposer.def( "__deepcopy__", &__copy__);
        MoleculeBeading_exposer.def( "clone", &__copy__);
        MoleculeBeading_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::MoleculeBeading >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MoleculeBeading_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::MoleculeBeading >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MoleculeBeading_exposer.def( "__str__", &__str__< ::SireMol::MoleculeBeading > );
        MoleculeBeading_exposer.def( "__repr__", &__str__< ::SireMol::MoleculeBeading > );
    }

}
