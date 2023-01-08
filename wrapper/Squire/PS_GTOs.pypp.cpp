// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "PS_GTOs.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include "SireMaths/boys.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "pgto.h"

#include "pointcharge.h"

#include "pointdipole.h"

#include "pgto.h"

Squire::PS_GTOs __copy__(const Squire::PS_GTOs &other){ return Squire::PS_GTOs(other); }

const char* pvt_get_name(const Squire::PS_GTOs&){ return "Squire::PS_GTOs";}

#include "Helpers/release_gil_policy.hpp"

void register_PS_GTOs_class(){

    { //::Squire::PS_GTOs
        typedef bp::class_< Squire::PS_GTOs > PS_GTOs_exposer_t;
        PS_GTOs_exposer_t PS_GTOs_exposer = PS_GTOs_exposer_t( "PS_GTOs", "This class is used to calculate integrals involving PS orbital\npairs\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope PS_GTOs_scope( PS_GTOs_exposer );
        PS_GTOs_exposer.def( bp::init< QVector< Squire::S_GTO > const &, QVector< SireMaths::Vector > const &, QVector< Squire::P_GTO > const &, QVector< SireMaths::Vector > const & >(( bp::arg("s_gtos"), bp::arg("s_centers"), bp::arg("p_gtos"), bp::arg("p_centers") ), "Construct for the set of s orbitals s_gtos located at centers\ns_centers, and p orbitals p_gtos located at centers\np_centers") );
        PS_GTOs_exposer.def( bp::init< QVector< Squire::P_GTO > const &, QVector< SireMaths::Vector > const &, QVector< Squire::S_GTO > const &, QVector< SireMaths::Vector > const & >(( bp::arg("p_gtos"), bp::arg("p_centers"), bp::arg("s_gtos"), bp::arg("s_centers") ), "Construct for the set of s orbitals s_gtos located at centers\ns_centers, and p orbitals p_gtos located at centers\np_centers") );
        PS_GTOs_exposer.def( bp::init< Squire::PS_GTOs const & >(( bp::arg("other") ), "Copy constructor") );
        { //::Squire::PS_GTOs::kinetic_integral
        
            typedef ::SireMaths::NMatrix ( ::Squire::PS_GTOs::*kinetic_integral_function_type)(  ) const;
            kinetic_integral_function_type kinetic_integral_function_value( &::Squire::PS_GTOs::kinetic_integral );
            
            PS_GTOs_exposer.def( 
                "kinetic_integral"
                , kinetic_integral_function_value
                , bp::release_gil_policy()
                , "Return the kinetic energy integrals for each pair of orbitals" );
        
        }
        PS_GTOs_exposer.def( bp::self != bp::self );
        { //::Squire::PS_GTOs::operator=
        
            typedef ::Squire::PS_GTOs & ( ::Squire::PS_GTOs::*assign_function_type)( ::Squire::PS_GTOs const & ) ;
            assign_function_type assign_function_value( &::Squire::PS_GTOs::operator= );
            
            PS_GTOs_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PS_GTOs_exposer.def( bp::self == bp::self );
        { //::Squire::PS_GTOs::overlap_integral
        
            typedef ::SireMaths::NMatrix ( ::Squire::PS_GTOs::*overlap_integral_function_type)(  ) const;
            overlap_integral_function_type overlap_integral_function_value( &::Squire::PS_GTOs::overlap_integral );
            
            PS_GTOs_exposer.def( 
                "overlap_integral"
                , overlap_integral_function_value
                , bp::release_gil_policy()
                , "Return the overlap integrals for each pair of orbitals" );
        
        }
        { //::Squire::PS_GTOs::potential_integral
        
            typedef ::SireMaths::NMatrix ( ::Squire::PS_GTOs::*potential_integral_function_type)( ::QVector< Squire::PointCharge > const & ) const;
            potential_integral_function_type potential_integral_function_value( &::Squire::PS_GTOs::potential_integral );
            
            PS_GTOs_exposer.def( 
                "potential_integral"
                , potential_integral_function_value
                , ( bp::arg("C") )
                , bp::release_gil_policy()
                , "Return the potential energy integral for each pair of orbitals\ninteracting with the point charges in C" );
        
        }
        { //::Squire::PS_GTOs::potential_integral
        
            typedef ::SireMaths::NMatrix ( ::Squire::PS_GTOs::*potential_integral_function_type)( ::QVector< Squire::PointCharge > const &,int ) const;
            potential_integral_function_type potential_integral_function_value( &::Squire::PS_GTOs::potential_integral );
            
            PS_GTOs_exposer.def( 
                "potential_integral"
                , potential_integral_function_value
                , ( bp::arg("C"), bp::arg("m") )
                , bp::release_gil_policy()
                , "Return the mth auxillary potential energy integral for each pair of orbitals\ninteracting with the point charges in C" );
        
        }
        PS_GTOs_exposer.def( "__copy__", &__copy__);
        PS_GTOs_exposer.def( "__deepcopy__", &__copy__);
        PS_GTOs_exposer.def( "clone", &__copy__);
        PS_GTOs_exposer.def( "__str__", &pvt_get_name);
        PS_GTOs_exposer.def( "__repr__", &pvt_get_name);
    }

}
