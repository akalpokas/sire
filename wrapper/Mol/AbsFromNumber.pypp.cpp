// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AbsFromNumber.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atomselection.h"

#include "editor.hpp"

#include "evaluator.h"

#include "moleculedata.h"

#include "moleculeview.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "weightfunction.h"

#include "weightfunction.h"

SireMol::AbsFromNumber __copy__(const SireMol::AbsFromNumber &other){ return SireMol::AbsFromNumber(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AbsFromNumber_class(){

    { //::SireMol::AbsFromNumber
        typedef bp::class_< SireMol::AbsFromNumber, bp::bases< SireMol::WeightFunction, SireBase::Property > > AbsFromNumber_exposer_t;
        AbsFromNumber_exposer_t AbsFromNumber_exposer = AbsFromNumber_exposer_t( "AbsFromNumber", "This class calculates the weight by assigning all of the weight to\nthe group with the largest number of atoms.\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope AbsFromNumber_scope( AbsFromNumber_exposer );
        AbsFromNumber_exposer.def( bp::init< SireMol::AbsFromNumber const & >(( bp::arg("other") ), "") );
        AbsFromNumber_exposer.def( bp::self != bp::self );
        { //::SireMol::AbsFromNumber::operator()
        
            typedef double ( ::SireMol::AbsFromNumber::*__call___function_type)( ::SireMol::MoleculeData const &,::SireMol::AtomSelection const &,::SireMol::AtomSelection const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::AbsFromNumber::operator() );
            
            AbsFromNumber_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("moldata"), bp::arg("group0"), bp::arg("group1"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::AbsFromNumber::operator()
        
            typedef double ( ::SireMol::AbsFromNumber::*__call___function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::AbsFromNumber::operator() );
            
            AbsFromNumber_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("view0"), bp::arg("map0"), bp::arg("view1"), bp::arg("map1") )
                , "" );
        
        }
        { //::SireMol::AbsFromNumber::operator=
        
            typedef ::SireMol::AbsFromNumber & ( ::SireMol::AbsFromNumber::*assign_function_type)( ::SireMol::AbsFromNumber const & ) ;
            assign_function_type assign_function_value( &::SireMol::AbsFromNumber::operator= );
            
            AbsFromNumber_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("arg0") )
                , bp::return_self< >()
                , "" );
        
        }
        AbsFromNumber_exposer.def( bp::self == bp::self );
        { //::SireMol::AbsFromNumber::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AbsFromNumber::typeName );
            
            AbsFromNumber_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        AbsFromNumber_exposer.staticmethod( "typeName" );
        AbsFromNumber_exposer.def( "__copy__", &__copy__);
        AbsFromNumber_exposer.def( "__deepcopy__", &__copy__);
        AbsFromNumber_exposer.def( "clone", &__copy__);
        AbsFromNumber_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AbsFromNumber >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AbsFromNumber_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AbsFromNumber >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AbsFromNumber_exposer.def_pickle(sire_pickle_suite< ::SireMol::AbsFromNumber >());
        AbsFromNumber_exposer.def( "__str__", &__str__< ::SireMol::AbsFromNumber > );
        AbsFromNumber_exposer.def( "__repr__", &__str__< ::SireMol::AbsFromNumber > );
    }

}
