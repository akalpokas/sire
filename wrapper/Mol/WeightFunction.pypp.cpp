// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "WeightFunction.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atomselection.h"

#include "core.h"

#include "editor.hpp"

#include "evaluator.h"

#include "moleculedata.h"

#include "moleculeview.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "weightfunction.h"

#include "weightfunction.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_WeightFunction_class(){

    { //::SireMol::WeightFunction
        typedef bp::class_< SireMol::WeightFunction, bp::bases< SireBase::Property >, boost::noncopyable > WeightFunction_exposer_t;
        WeightFunction_exposer_t WeightFunction_exposer = WeightFunction_exposer_t( "WeightFunction", "This is the base class of all weight functions. A weight function\nis a simple function that takes two groups in a molecule, and\nreturns the relative weight of those two groups (0 == 100% group A,\n0.5 == 50% group A, 50% group B, 1 == 100% groupB)\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope WeightFunction_scope( WeightFunction_exposer );
        { //::SireMol::WeightFunction::null
        
            typedef ::SireMol::AbsFromNumber const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::WeightFunction::null );
            
            WeightFunction_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMol::WeightFunction::operator()
        
            typedef double ( ::SireMol::WeightFunction::*__call___function_type)( ::SireMol::MoleculeData const &,::SireMol::AtomSelection const &,::SireMol::AtomSelection const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::WeightFunction::operator() );
            
            WeightFunction_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("moldata"), bp::arg("group0"), bp::arg("group1"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the relative weight of group0 and group1 in the\nmolecule whose data is in moldata, using the supplied\nmap to find the required properties\n\nThrow: SireError::incompatible_error\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::WeightFunction::operator()
        
            typedef double ( ::SireMol::WeightFunction::*__call___function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::WeightFunction::operator() );
            
            WeightFunction_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("view0"), bp::arg("map0"), bp::arg("view1"), bp::arg("map1") )
                , "Return the relative weight of two molecule views (view0\nand view1), using map0 to find the required properties\nfrom view0, and map1 to find the required properties from\nview1.\n\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::WeightFunction::operator()
        
            typedef double ( ::SireMol::WeightFunction::*__call___function_type)( ::SireMol::MoleculeView const &,::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::WeightFunction::operator() );
            
            WeightFunction_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("view0"), bp::arg("view1"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the relative weight of the two molecule views\n(view0 and view1) using the supplied map to find the\nrequired properties from both views\n\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::WeightFunction::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::WeightFunction::typeName );
            
            WeightFunction_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        WeightFunction_exposer.staticmethod( "null" );
        WeightFunction_exposer.staticmethod( "typeName" );
        WeightFunction_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::WeightFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        WeightFunction_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::WeightFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        WeightFunction_exposer.def_pickle(sire_pickle_suite< ::SireMol::WeightFunction >());
        WeightFunction_exposer.def( "__str__", &__str__< ::SireMol::WeightFunction > );
        WeightFunction_exposer.def( "__repr__", &__str__< ::SireMol::WeightFunction > );
    }

}
