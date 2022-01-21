// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PMF.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/registeralternativename.h"

#include "SireStream/shareddatastream.h"

#include "fep.h"

#include "tostring.h"

#include "fep.h"

SireAnalysis::PMF __copy__(const SireAnalysis::PMF &other){ return SireAnalysis::PMF(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_PMF_class(){

    { //::SireAnalysis::PMF
        typedef bp::class_< SireAnalysis::PMF, bp::bases< SireBase::Property > > PMF_exposer_t;
        PMF_exposer_t PMF_exposer = PMF_exposer_t( "PMF", "This class contains the complete potential of mean force\n\nAuthor: Christopher Woods\n", bp::init< >("Construct an empty PMF") );
        bp::scope PMF_scope( PMF_exposer );
        PMF_exposer.def( bp::init< QVector< SireAnalysis::DataPoint > const & >(( bp::arg("values") ), "Construct a from the passed data points") );
        PMF_exposer.def( bp::init< SireAnalysis::PMF const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireAnalysis::PMF::deltaG
        
            typedef double ( ::SireAnalysis::PMF::*deltaG_function_type)(  ) const;
            deltaG_function_type deltaG_function_value( &::SireAnalysis::PMF::deltaG );
            
            PMF_exposer.def( 
                "deltaG"
                , deltaG_function_value
                , "Return the total free energy change along the PMF (difference in\nfree energy of the end-points)" );
        
        }
        { //::SireAnalysis::PMF::error
        
            typedef double ( ::SireAnalysis::PMF::*error_function_type)(  ) const;
            error_function_type error_function_value( &::SireAnalysis::PMF::error );
            
            PMF_exposer.def( 
                "error"
                , error_function_value
                , "Return the error on the total free energy calculation" );
        
        }
        { //::SireAnalysis::PMF::isEmpty
        
            typedef bool ( ::SireAnalysis::PMF::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireAnalysis::PMF::isEmpty );
            
            PMF_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , "Return whether or not this PMF is empty (has not values)" );
        
        }
        PMF_exposer.def( bp::self != bp::self );
        { //::SireAnalysis::PMF::operator=
        
            typedef ::SireAnalysis::PMF & ( ::SireAnalysis::PMF::*assign_function_type)( ::SireAnalysis::PMF const & ) ;
            assign_function_type assign_function_value( &::SireAnalysis::PMF::operator= );
            
            PMF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PMF_exposer.def( bp::self == bp::self );
        { //::SireAnalysis::PMF::rangeMax
        
            typedef double ( ::SireAnalysis::PMF::*rangeMax_function_type)(  ) const;
            rangeMax_function_type rangeMax_function_value( &::SireAnalysis::PMF::rangeMax );
            
            PMF_exposer.def( 
                "rangeMax"
                , rangeMax_function_value
                , "Return the maximum x-value of the PMF" );
        
        }
        { //::SireAnalysis::PMF::rangeMin
        
            typedef double ( ::SireAnalysis::PMF::*rangeMin_function_type)(  ) const;
            rangeMin_function_type rangeMin_function_value( &::SireAnalysis::PMF::rangeMin );
            
            PMF_exposer.def( 
                "rangeMin"
                , rangeMin_function_value
                , "Return the minimum x-value of the PMF" );
        
        }
        { //::SireAnalysis::PMF::toString
        
            typedef ::QString ( ::SireAnalysis::PMF::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireAnalysis::PMF::toString );
            
            PMF_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireAnalysis::PMF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireAnalysis::PMF::typeName );
            
            PMF_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireAnalysis::PMF::values
        
            typedef ::QVector< SireAnalysis::DataPoint > ( ::SireAnalysis::PMF::*values_function_type)(  ) const;
            values_function_type values_function_value( &::SireAnalysis::PMF::values );
            
            PMF_exposer.def( 
                "values"
                , values_function_value
                , "Return the raw data for the PMF" );
        
        }
        { //::SireAnalysis::PMF::what
        
            typedef char const * ( ::SireAnalysis::PMF::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireAnalysis::PMF::what );
            
            PMF_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        PMF_exposer.staticmethod( "typeName" );
        PMF_exposer.def( "__copy__", &__copy__);
        PMF_exposer.def( "__deepcopy__", &__copy__);
        PMF_exposer.def( "clone", &__copy__);
        PMF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireAnalysis::PMF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PMF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireAnalysis::PMF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PMF_exposer.def_pickle(sire_pickle_suite< ::SireAnalysis::PMF >());
        PMF_exposer.def( "__str__", &__str__< ::SireAnalysis::PMF > );
        PMF_exposer.def( "__repr__", &__str__< ::SireAnalysis::PMF > );
    }

}
