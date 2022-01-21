// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CovalentBondHunter.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/coordgroup.h"

#include "atom.h"

#include "atomcoords.h"

#include "atomelements.h"

#include "atomselection.h"

#include "bondhunter.h"

#include "connectivity.h"

#include "molecule.h"

#include "moleculedata.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "mover.hpp"

#include "selector.hpp"

#include <QDebug>

#include <QMutex>

#include "bondhunter.h"

SireMol::CovalentBondHunter __copy__(const SireMol::CovalentBondHunter &other){ return SireMol::CovalentBondHunter(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CovalentBondHunter_class(){

    { //::SireMol::CovalentBondHunter
        typedef bp::class_< SireMol::CovalentBondHunter, bp::bases< SireMol::BondHunter, SireBase::Property > > CovalentBondHunter_exposer_t;
        CovalentBondHunter_exposer_t CovalentBondHunter_exposer = CovalentBondHunter_exposer_t( "CovalentBondHunter", "This is a bond hunter that finds bonded atoms by comparing\nthe distance between the atoms to the sum of their covalent radii.\nIf the distance is less than the sum of the covalent radii then the atoms\nare bonded.\n\nAuthor: Christopher Woods\n", bp::init< bp::optional< double > >(( bp::arg("tolerance")=1.1000000000000001 ), "Construct a CovalentBondHunter with a specified tolerance") );
        bp::scope CovalentBondHunter_scope( CovalentBondHunter_exposer );
        CovalentBondHunter_exposer.def( bp::init< SireMol::CovalentBondHunter const & >(( bp::arg("other") ), "Copy constructor") );
        CovalentBondHunter_exposer.def( bp::self != bp::self );
        { //::SireMol::CovalentBondHunter::operator()
        
            typedef ::SireMol::Connectivity ( ::SireMol::CovalentBondHunter::*__call___function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::CovalentBondHunter::operator() );
            
            CovalentBondHunter_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::CovalentBondHunter::operator=
        
            typedef ::SireMol::CovalentBondHunter & ( ::SireMol::CovalentBondHunter::*assign_function_type)( ::SireMol::CovalentBondHunter const & ) ;
            assign_function_type assign_function_value( &::SireMol::CovalentBondHunter::operator= );
            
            CovalentBondHunter_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CovalentBondHunter_exposer.def( bp::self == bp::self );
        { //::SireMol::CovalentBondHunter::parameters
        
            typedef ::SireMol::CovalentBondHunterParameters const & ( *parameters_function_type )(  );
            parameters_function_type parameters_function_value( &::SireMol::CovalentBondHunter::parameters );
            
            CovalentBondHunter_exposer.def( 
                "parameters"
                , parameters_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::CovalentBondHunter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CovalentBondHunter::typeName );
            
            CovalentBondHunter_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        CovalentBondHunter_exposer.staticmethod( "parameters" );
        CovalentBondHunter_exposer.staticmethod( "typeName" );
        CovalentBondHunter_exposer.def( "__copy__", &__copy__);
        CovalentBondHunter_exposer.def( "__deepcopy__", &__copy__);
        CovalentBondHunter_exposer.def( "clone", &__copy__);
        CovalentBondHunter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::CovalentBondHunter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CovalentBondHunter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::CovalentBondHunter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CovalentBondHunter_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::CovalentBondHunter > );
        CovalentBondHunter_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::CovalentBondHunter > );
        CovalentBondHunter_exposer.def( "__str__", &__str__< ::SireMol::CovalentBondHunter > );
        CovalentBondHunter_exposer.def( "__repr__", &__str__< ::SireMol::CovalentBondHunter > );
    }

}
