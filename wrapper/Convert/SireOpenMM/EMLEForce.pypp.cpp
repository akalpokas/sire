// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "EMLEForce.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireVol/triclinicbox.h"

#include "emle.h"

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireVol/triclinicbox.h"

#include "emle.h"

SireOpenMM::EMLEForce __copy__(const SireOpenMM::EMLEForce &other){ return SireOpenMM::EMLEForce(other); }

const char* pvt_get_name(const SireOpenMM::EMLEForce&){ return "SireOpenMM::EMLEForce";}

#include "Helpers/release_gil_policy.hpp"

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireVol/triclinicbox.h"

#include "emle.h"

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireVol/triclinicbox.h"

#include "emle.h"

const char* pvt_get_name(const SireOpenMM::EMLEForce&){ return "SireOpenMM::EMLEForce";}

#include "Helpers/release_gil_policy.hpp"

void register_EMLEForce_class(){

    { //::SireOpenMM::EMLEForce
        typedef bp::class_< SireOpenMM::EMLEForce, bp::bases< SireOpenMM::QMForce > > EMLEForce_exposer_t;
        EMLEForce_exposer_t EMLEForce_exposer = EMLEForce_exposer_t( "EMLEForce", "", bp::init< >("Default constructor.") );
        bp::scope EMLEForce_scope( EMLEForce_exposer );
        EMLEForce_exposer.def( bp::init< SireOpenMM::EMLECallback, SireUnits::Dimension::Length, int, double, QVector< int >, QMap< int, int >, QMap< int, QVector< int > >, QMap< int, double >, QVector< int >, QVector< int >, QVector< double > >(( bp::arg("callback"), bp::arg("cutoff"), bp::arg("neighbour_list_frequency"), bp::arg("lambda"), bp::arg("atoms"), bp::arg("mm1_to_qm"), bp::arg("mm1_to_mm2"), bp::arg("bond_scale_factors"), bp::arg("mm2_atoms"), bp::arg("numbers"), bp::arg("charges") ), "Constructor.\nPar:am callback\nThe EMLECallback object.\n\nPar:am cutoff\nThe ML cutoff distance.\n\nPar:am neighbour_list_frequency\nThe frequency at which the neighbour list is updated. (Number of steps.)\nIf zero, then no neighbour list is used.\n\nPar:am lambda\nThe lambda weighting factor. This can be used to interpolate between\npotentials for end-state correction calculations.\n\nPar:am atoms\nA vector of atom indices for the QM region.\n\nPar:am mm1_to_qm\nA dictionary mapping link atom (MM1) indices to the QM atoms to\nwhich they are bonded.\n\nPar:am mm1_to_mm2\nA dictionary of link atoms indices (MM1) to a list of the MM\natoms to which they are bonded (MM2).\n\nPar:am bond_scale_factors\nA dictionary of link atom indices (MM1) to a list of the bond\nlength scale factors between the QM and MM1 atoms. The scale\nfactors are the ratio of the equilibrium bond lengths for the\nQM-L (QM-link) atom and QM-MM1 atom, i.e. R0(QM-L)  R0(QM-MM1),\ntaken from the MM force field parameters for the molecule.\n\nPar:am mm2_atoms\nA vector of MM2 atom indices.\n\nPar:am numbers\nA vector of atomic charges for all atoms in the system.\n\nPar:am charges\nA vector of atomic charges for all atoms in the system.\n") );
        EMLEForce_exposer.def( bp::init< SireOpenMM::EMLEForce const & >(( bp::arg("other") ), "Copy constructor.") );
        { //::SireOpenMM::EMLEForce::call
        
            typedef ::boost::tuples::tuple< double, QVector< QVector< double > >, QVector< QVector< double > >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireOpenMM::EMLEForce::*call_function_type)( ::QVector< int >,::QVector< double >,::QVector< QVector< double > >,::QVector< QVector< double > > ) const;
            call_function_type call_function_value( &::SireOpenMM::EMLEForce::call );
            
            EMLEForce_exposer.def( 
                "call"
                , call_function_value
                , ( bp::arg("numbers_qm"), bp::arg("charges_mm"), bp::arg("xyz_qm"), bp::arg("xyz_mm") )
                , bp::release_gil_policy()
                , "Call the callback function.\nPar:am numbers_qm\nA vector of atomic numbers for the atoms in the ML region.\n\nPar:am charges_mm\nA vector of the charges on the MM atoms in mod electron charge.\n\nPar:am xyz_qm\nA vector of positions for the atoms in the ML region in Angstrom.\n\nPar:am xyz_mm\nA vector of positions for the atoms in the MM region in Angstrom.\n\nReturn:s\nA tuple containing:\n- The energy in kJmol.\n- A vector of forces for the QM atoms in kJmolnm.\n- A vector of forces for the MM atoms in kJmolnm.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getAtoms
        
            typedef ::QVector< int > ( ::SireOpenMM::EMLEForce::*getAtoms_function_type)(  ) const;
            getAtoms_function_type getAtoms_function_value( &::SireOpenMM::EMLEForce::getAtoms );
            
            EMLEForce_exposer.def( 
                "getAtoms"
                , getAtoms_function_value
                , bp::release_gil_policy()
                , "Get the indices of the atoms in the QM region.\nReturn:s\nA vector of atom indices for the QM region.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getCallback
        
            typedef ::SireOpenMM::EMLECallback ( ::SireOpenMM::EMLEForce::*getCallback_function_type)(  ) const;
            getCallback_function_type getCallback_function_value( &::SireOpenMM::EMLEForce::getCallback );
            
            EMLEForce_exposer.def( 
                "getCallback"
                , getCallback_function_value
                , bp::release_gil_policy()
                , "Get the callback object.\nReturn:s\nA Python object that contains the callback function.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getCharges
        
            typedef ::QVector< double > ( ::SireOpenMM::EMLEForce::*getCharges_function_type)(  ) const;
            getCharges_function_type getCharges_function_value( &::SireOpenMM::EMLEForce::getCharges );
            
            EMLEForce_exposer.def( 
                "getCharges"
                , getCharges_function_value
                , bp::release_gil_policy()
                , "Get the atomic charges of all atoms in the system.\nReturn:s\nA vector of atomic charges for all atoms in the system.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getCutoff
        
            typedef ::SireUnits::Dimension::Length ( ::SireOpenMM::EMLEForce::*getCutoff_function_type)(  ) const;
            getCutoff_function_type getCutoff_function_value( &::SireOpenMM::EMLEForce::getCutoff );
            
            EMLEForce_exposer.def( 
                "getCutoff"
                , getCutoff_function_value
                , bp::release_gil_policy()
                , "Get the QM cutoff distance.\nReturn:s\nThe QM cutoff distance.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getLambda
        
            typedef double ( ::SireOpenMM::EMLEForce::*getLambda_function_type)(  ) const;
            getLambda_function_type getLambda_function_value( &::SireOpenMM::EMLEForce::getLambda );
            
            EMLEForce_exposer.def( 
                "getLambda"
                , getLambda_function_value
                , bp::release_gil_policy()
                , "Get the lambda weighting factor.\nReturn:s\nThe lambda weighting factor.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getLinkAtoms
        
            typedef ::boost::tuples::tuple< QMap< int, int >, QMap< int, QVector< int > >, QMap< int, double >, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireOpenMM::EMLEForce::*getLinkAtoms_function_type)(  ) const;
            getLinkAtoms_function_type getLinkAtoms_function_value( &::SireOpenMM::EMLEForce::getLinkAtoms );
            
            EMLEForce_exposer.def( 
                "getLinkAtoms"
                , getLinkAtoms_function_value
                , bp::release_gil_policy()
                , "Get the link atoms associated with each QM atom.\nReturn:s\nA tuple containing:\n\nmm1_to_qm\nA dictionary mapping link atom (MM1) indices to the QM atoms to\nwhich they are bonded.\n\nmm1_to_mm2\nA dictionary of link atoms indices (MM1) to a list of the MM\natoms to which they are bonded (MM2).\n\nbond_scale_factors\nA dictionary of link atom indices (MM1) to a list of the bond\nlength scale factors between the QM and MM1 atoms. The scale\nfactors are the ratio of the equilibrium bond lengths for the\nQM-L (QM-link) atom and QM-MM1 atom, i.e. R0(QM-L)  R0(QM-MM1),\ntaken from the MM force field parameters for the molecule.\n\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getMM2Atoms
        
            typedef ::QVector< int > ( ::SireOpenMM::EMLEForce::*getMM2Atoms_function_type)(  ) const;
            getMM2Atoms_function_type getMM2Atoms_function_value( &::SireOpenMM::EMLEForce::getMM2Atoms );
            
            EMLEForce_exposer.def( 
                "getMM2Atoms"
                , getMM2Atoms_function_value
                , bp::release_gil_policy()
                , "Get the vector of MM2 atoms.\nReturn:s\nA vector of MM2 atom indices.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getNeighbourListFrequency
        
            typedef int ( ::SireOpenMM::EMLEForce::*getNeighbourListFrequency_function_type)(  ) const;
            getNeighbourListFrequency_function_type getNeighbourListFrequency_function_value( &::SireOpenMM::EMLEForce::getNeighbourListFrequency );
            
            EMLEForce_exposer.def( 
                "getNeighbourListFrequency"
                , getNeighbourListFrequency_function_value
                , bp::release_gil_policy()
                , "Get the neighbour list frequency.\nReturn:s\nThe neighbour list frequency.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::getNumbers
        
            typedef ::QVector< int > ( ::SireOpenMM::EMLEForce::*getNumbers_function_type)(  ) const;
            getNumbers_function_type getNumbers_function_value( &::SireOpenMM::EMLEForce::getNumbers );
            
            EMLEForce_exposer.def( 
                "getNumbers"
                , getNumbers_function_value
                , bp::release_gil_policy()
                , "Get the atomic numbers for the atoms in the QM region.\nReturn:s\nA vector of atomic numbers for the atoms in the QM region.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::operator=
        
            typedef ::SireOpenMM::EMLEForce & ( ::SireOpenMM::EMLEForce::*assign_function_type)( ::SireOpenMM::EMLEForce const & ) ;
            assign_function_type assign_function_value( &::SireOpenMM::EMLEForce::operator= );
            
            EMLEForce_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "Assignment operator." );
        
        }
        { //::SireOpenMM::EMLEForce::setLambda
        
            typedef void ( ::SireOpenMM::EMLEForce::*setLambda_function_type)( double ) ;
            setLambda_function_type setLambda_function_value( &::SireOpenMM::EMLEForce::setLambda );
            
            EMLEForce_exposer.def( 
                "setLambda"
                , setLambda_function_value
                , ( bp::arg("lambda") )
                , bp::release_gil_policy()
                , "Set the lambda weighting factor\nPar:am lambda\nThe lambda weighting factor.\n" );
        
        }
        { //::SireOpenMM::EMLEForce::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireOpenMM::EMLEForce::typeName );
            
            EMLEForce_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "Return the C++ name for this class." );
        
        }
        { //::SireOpenMM::EMLEForce::what
        
            typedef char const * ( ::SireOpenMM::EMLEForce::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireOpenMM::EMLEForce::what );
            
            EMLEForce_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "Return the C++ name for this class." );
        
        }
        EMLEForce_exposer.staticmethod( "typeName" );
        EMLEForce_exposer.def( "__copy__", &__copy__);
        EMLEForce_exposer.def( "__deepcopy__", &__copy__);
        EMLEForce_exposer.def( "clone", &__copy__);
        EMLEForce_exposer.def( "__str__", &pvt_get_name);
        EMLEForce_exposer.def( "__repr__", &pvt_get_name);
        EMLEForce_exposer.def( "__str__", &pvt_get_name);
        EMLEForce_exposer.def( "__repr__", &pvt_get_name);
    }

}
