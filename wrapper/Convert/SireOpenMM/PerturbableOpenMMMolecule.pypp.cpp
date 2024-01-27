// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "PerturbableOpenMMMolecule.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/parallel.h"

#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireMM/amberparams.h"

#include "SireMM/atomljs.h"

#include "SireMM/selectorbond.h"

#include "SireMM/twoatomfunctions.h"

#include "SireMaths/vector.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/atomproperty.hpp"

#include "SireMol/atomvelocities.h"

#include "SireMol/bondid.h"

#include "SireMol/bondorder.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/moleditor.h"

#include "SireUnits/units.h"

#include "openmmmolecule.h"

#include "tostring.h"

#include <QDebug>

#include <QReadWriteLock>

#include <QSet>

#include "SireBase/parallel.h"

#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireMM/amberparams.h"

#include "SireMM/atomljs.h"

#include "SireMM/selectorbond.h"

#include "SireMM/twoatomfunctions.h"

#include "SireMaths/vector.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/atomproperty.hpp"

#include "SireMol/atomvelocities.h"

#include "SireMol/bondid.h"

#include "SireMol/bondorder.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/moleditor.h"

#include "SireUnits/units.h"

#include "openmmmolecule.h"

#include "tostring.h"

#include <QDebug>

#include <QReadWriteLock>

#include <QSet>

SireOpenMM::PerturbableOpenMMMolecule __copy__(const SireOpenMM::PerturbableOpenMMMolecule &other){ return SireOpenMM::PerturbableOpenMMMolecule(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "SireBase/parallel.h"

#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireMM/amberparams.h"

#include "SireMM/atomljs.h"

#include "SireMM/selectorbond.h"

#include "SireMM/twoatomfunctions.h"

#include "SireMaths/vector.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/atomproperty.hpp"

#include "SireMol/atomvelocities.h"

#include "SireMol/bondid.h"

#include "SireMol/bondorder.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/moleditor.h"

#include "SireUnits/units.h"

#include "openmmmolecule.h"

#include "tostring.h"

#include <QDebug>

#include <QReadWriteLock>

#include <QSet>

#include "SireBase/parallel.h"

#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireMM/amberparams.h"

#include "SireMM/atomljs.h"

#include "SireMM/selectorbond.h"

#include "SireMM/twoatomfunctions.h"

#include "SireMaths/vector.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/atomproperty.hpp"

#include "SireMol/atomvelocities.h"

#include "SireMol/bondid.h"

#include "SireMol/bondorder.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/moleditor.h"

#include "SireUnits/units.h"

#include "openmmmolecule.h"

#include "tostring.h"

#include <QDebug>

#include <QReadWriteLock>

#include <QSet>

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_PerturbableOpenMMMolecule_class(){

    { //::SireOpenMM::PerturbableOpenMMMolecule
        typedef bp::class_< SireOpenMM::PerturbableOpenMMMolecule > PerturbableOpenMMMolecule_exposer_t;
        PerturbableOpenMMMolecule_exposer_t PerturbableOpenMMMolecule_exposer = PerturbableOpenMMMolecule_exposer_t( "PerturbableOpenMMMolecule", "This class holds all of the information of an OpenMM molecule\nthat can be perturbed using a LambdaSchedule. The data is held\nin easy-to-access arrays, with guarantees that the arrays are\ncompatible and the data is aligned.\n", bp::init< >("Null constructor") );
        bp::scope PerturbableOpenMMMolecule_scope( PerturbableOpenMMMolecule_exposer );
        PerturbableOpenMMMolecule_exposer.def( bp::init< SireOpenMM::OpenMMMolecule const & >(( bp::arg("mol") ), "Construct from the passed OpenMMMolecule") );
        PerturbableOpenMMMolecule_exposer.def( bp::init< SireOpenMM::PerturbableOpenMMMolecule const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireOpenMM::PerturbableOpenMMMolecule::getAlphas0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getAlphas0_function_type)(  ) const;
            getAlphas0_function_type getAlphas0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getAlphas0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getAlphas0"
                , getAlphas0_function_value
                , bp::release_gil_policy()
                , "Return the alpha parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getAlphas1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getAlphas1_function_type)(  ) const;
            getAlphas1_function_type getAlphas1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getAlphas1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getAlphas1"
                , getAlphas1_function_value
                , bp::release_gil_policy()
                , "Return the alpha parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getAngleKs0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getAngleKs0_function_type)(  ) const;
            getAngleKs0_function_type getAngleKs0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getAngleKs0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getAngleKs0"
                , getAngleKs0_function_value
                , bp::release_gil_policy()
                , "Return the angle k parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getAngleKs1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getAngleKs1_function_type)(  ) const;
            getAngleKs1_function_type getAngleKs1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getAngleKs1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getAngleKs1"
                , getAngleKs1_function_value
                , bp::release_gil_policy()
                , "Return the angle k parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getAngleSizes0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getAngleSizes0_function_type)(  ) const;
            getAngleSizes0_function_type getAngleSizes0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getAngleSizes0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getAngleSizes0"
                , getAngleSizes0_function_value
                , bp::release_gil_policy()
                , "Return the angle size parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getAngleSizes1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getAngleSizes1_function_type)(  ) const;
            getAngleSizes1_function_type getAngleSizes1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getAngleSizes1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getAngleSizes1"
                , getAngleSizes1_function_value
                , bp::release_gil_policy()
                , "Return the angle size parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getBondKs0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getBondKs0_function_type)(  ) const;
            getBondKs0_function_type getBondKs0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getBondKs0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getBondKs0"
                , getBondKs0_function_value
                , bp::release_gil_policy()
                , "Return the bond k parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getBondKs1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getBondKs1_function_type)(  ) const;
            getBondKs1_function_type getBondKs1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getBondKs1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getBondKs1"
                , getBondKs1_function_value
                , bp::release_gil_policy()
                , "Return the bond k parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getBondLengths0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getBondLengths0_function_type)(  ) const;
            getBondLengths0_function_type getBondLengths0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getBondLengths0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getBondLengths0"
                , getBondLengths0_function_value
                , bp::release_gil_policy()
                , "Return the bond length parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getBondLengths1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getBondLengths1_function_type)(  ) const;
            getBondLengths1_function_type getBondLengths1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getBondLengths1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getBondLengths1"
                , getBondLengths1_function_value
                , bp::release_gil_policy()
                , "Return the bond length parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getChargeScales0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getChargeScales0_function_type)(  ) const;
            getChargeScales0_function_type getChargeScales0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getChargeScales0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getChargeScales0"
                , getChargeScales0_function_value
                , bp::release_gil_policy()
                , "Return the coulomb intramolecular scale factors of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getChargeScales1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getChargeScales1_function_type)(  ) const;
            getChargeScales1_function_type getChargeScales1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getChargeScales1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getChargeScales1"
                , getChargeScales1_function_value
                , bp::release_gil_policy()
                , "Return the coulomb intramolecular scale factors of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getCharges0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getCharges0_function_type)(  ) const;
            getCharges0_function_type getCharges0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getCharges0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getCharges0"
                , getCharges0_function_value
                , bp::release_gil_policy()
                , "Return the atom charges of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getCharges1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getCharges1_function_type)(  ) const;
            getCharges1_function_type getCharges1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getCharges1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getCharges1"
                , getCharges1_function_value
                , bp::release_gil_policy()
                , "Return the atom charges of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getConstraintIndicies
        
            typedef ::QVector< int > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getConstraintIndicies_function_type)(  ) const;
            getConstraintIndicies_function_type getConstraintIndicies_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getConstraintIndicies );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getConstraintIndicies"
                , getConstraintIndicies_function_value
                , bp::release_gil_policy()
                , "Return the indicies of the perturbable constraints" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getEpsilons0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getEpsilons0_function_type)(  ) const;
            getEpsilons0_function_type getEpsilons0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getEpsilons0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getEpsilons0"
                , getEpsilons0_function_value
                , bp::release_gil_policy()
                , "Return the LJ epsilon parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getEpsilons1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getEpsilons1_function_type)(  ) const;
            getEpsilons1_function_type getEpsilons1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getEpsilons1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getEpsilons1"
                , getEpsilons1_function_value
                , bp::release_gil_policy()
                , "Return the LJ epsilon parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getExceptionAtoms
        
            typedef ::QVector< std::pair< int, int > > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getExceptionAtoms_function_type)(  ) const;
            getExceptionAtoms_function_type getExceptionAtoms_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getExceptionAtoms );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getExceptionAtoms"
                , getExceptionAtoms_function_value
                , bp::release_gil_policy()
                , "Return the indices of the atoms in the exceptions" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getExceptionIndicies
        
            typedef ::QVector< std::pair< int, int > > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getExceptionIndicies_function_type)( ::QString const & ) const;
            getExceptionIndicies_function_type getExceptionIndicies_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getExceptionIndicies );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getExceptionIndicies"
                , getExceptionIndicies_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Return the global indexes of the exceptions in the non-bonded and\n  ghost-14 forces\n" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getFromGhostIdxs
        
            typedef ::QSet< int > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getFromGhostIdxs_function_type)(  ) const;
            getFromGhostIdxs_function_type getFromGhostIdxs_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getFromGhostIdxs );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getFromGhostIdxs"
                , getFromGhostIdxs_function_value
                , bp::release_gil_policy()
                , "Return the indexes of the atoms that were ghosts in the\n  reference state\n" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getKappas0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getKappas0_function_type)(  ) const;
            getKappas0_function_type getKappas0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getKappas0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getKappas0"
                , getKappas0_function_value
                , bp::release_gil_policy()
                , "Return the kappa parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getKappas1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getKappas1_function_type)(  ) const;
            getKappas1_function_type getKappas1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getKappas1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getKappas1"
                , getKappas1_function_value
                , bp::release_gil_policy()
                , "Return the kappa parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getLJScales0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getLJScales0_function_type)(  ) const;
            getLJScales0_function_type getLJScales0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getLJScales0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getLJScales0"
                , getLJScales0_function_value
                , bp::release_gil_policy()
                , "Return the LJ intramolecular scale factors of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getLJScales1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getLJScales1_function_type)(  ) const;
            getLJScales1_function_type getLJScales1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getLJScales1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getLJScales1"
                , getLJScales1_function_value
                , bp::release_gil_policy()
                , "Return the LJ intramolecular scale factors of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getPerturbableConstraints
        
            typedef ::std::tuple< QVector< int >, QVector< double >, QVector< double > > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getPerturbableConstraints_function_type)(  ) const;
            getPerturbableConstraints_function_type getPerturbableConstraints_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getPerturbableConstraints );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getPerturbableConstraints"
                , getPerturbableConstraints_function_value
                , bp::release_gil_policy()
                , "Return three arrays containing the constraint indexes, and the\n  reference and perturbed values of the constraint lengths\n" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getSigmas0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getSigmas0_function_type)(  ) const;
            getSigmas0_function_type getSigmas0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getSigmas0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getSigmas0"
                , getSigmas0_function_value
                , bp::release_gil_policy()
                , "Return the LJ sigma parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getSigmas1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getSigmas1_function_type)(  ) const;
            getSigmas1_function_type getSigmas1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getSigmas1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getSigmas1"
                , getSigmas1_function_value
                , bp::release_gil_policy()
                , "Return the LJ sigma parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getToGhostIdxs
        
            typedef ::QSet< int > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getToGhostIdxs_function_type)(  ) const;
            getToGhostIdxs_function_type getToGhostIdxs_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getToGhostIdxs );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getToGhostIdxs"
                , getToGhostIdxs_function_value
                , bp::release_gil_policy()
                , "Return the indexes of the atoms that are to be ghosted in the\n  perturbed state\n" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getTorsionKs0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getTorsionKs0_function_type)(  ) const;
            getTorsionKs0_function_type getTorsionKs0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getTorsionKs0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getTorsionKs0"
                , getTorsionKs0_function_value
                , bp::release_gil_policy()
                , "Return the torsion k parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getTorsionKs1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getTorsionKs1_function_type)(  ) const;
            getTorsionKs1_function_type getTorsionKs1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getTorsionKs1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getTorsionKs1"
                , getTorsionKs1_function_value
                , bp::release_gil_policy()
                , "Return the torsion k parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getTorsionPeriodicities0
        
            typedef ::QVector< signed char > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getTorsionPeriodicities0_function_type)(  ) const;
            getTorsionPeriodicities0_function_type getTorsionPeriodicities0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getTorsionPeriodicities0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getTorsionPeriodicities0"
                , getTorsionPeriodicities0_function_value
                , bp::release_gil_policy()
                , "Return the torsion periodicity parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getTorsionPeriodicities1
        
            typedef ::QVector< signed char > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getTorsionPeriodicities1_function_type)(  ) const;
            getTorsionPeriodicities1_function_type getTorsionPeriodicities1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getTorsionPeriodicities1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getTorsionPeriodicities1"
                , getTorsionPeriodicities1_function_value
                , bp::release_gil_policy()
                , "Return the torsion periodicity parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getTorsionPhases0
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getTorsionPhases0_function_type)(  ) const;
            getTorsionPhases0_function_type getTorsionPhases0_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getTorsionPhases0 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getTorsionPhases0"
                , getTorsionPhases0_function_value
                , bp::release_gil_policy()
                , "Return the torsion phase parameters of the reference state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::getTorsionPhases1
        
            typedef ::QVector< double > ( ::SireOpenMM::PerturbableOpenMMMolecule::*getTorsionPhases1_function_type)(  ) const;
            getTorsionPhases1_function_type getTorsionPhases1_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::getTorsionPhases1 );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "getTorsionPhases1"
                , getTorsionPhases1_function_value
                , bp::release_gil_policy()
                , "Return the torsion phase parameters of the perturbed state" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::isGhostAtom
        
            typedef bool ( ::SireOpenMM::PerturbableOpenMMMolecule::*isGhostAtom_function_type)( int ) const;
            isGhostAtom_function_type isGhostAtom_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::isGhostAtom );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "isGhostAtom"
                , isGhostAtom_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Return true if the atom is a ghost atom in the\n  referenece or perturbed states" );
        
        }
        PerturbableOpenMMMolecule_exposer.def( bp::self != bp::self );
        { //::SireOpenMM::PerturbableOpenMMMolecule::operator=
        
            typedef ::SireOpenMM::PerturbableOpenMMMolecule & ( ::SireOpenMM::PerturbableOpenMMMolecule::*assign_function_type)( ::SireOpenMM::PerturbableOpenMMMolecule const & ) ;
            assign_function_type assign_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::operator= );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PerturbableOpenMMMolecule_exposer.def( bp::self == bp::self );
        { //::SireOpenMM::PerturbableOpenMMMolecule::setConstraintIndicies
        
            typedef void ( ::SireOpenMM::PerturbableOpenMMMolecule::*setConstraintIndicies_function_type)( ::QVector< int > const & ) ;
            setConstraintIndicies_function_type setConstraintIndicies_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::setConstraintIndicies );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "setConstraintIndicies"
                , setConstraintIndicies_function_value
                , ( bp::arg("constraint_idxs") )
                , bp::release_gil_policy()
                , "Set the indexes of perturbable constraints in the System" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::setExceptionIndicies
        
            typedef void ( ::SireOpenMM::PerturbableOpenMMMolecule::*setExceptionIndicies_function_type)( ::QString const &,::QVector< std::pair< int, int > > const & ) ;
            setExceptionIndicies_function_type setExceptionIndicies_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::setExceptionIndicies );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "setExceptionIndicies"
                , setExceptionIndicies_function_value
                , ( bp::arg("name"), bp::arg("exception_idxs") )
                , bp::release_gil_policy()
                , "Set the global indexes of the exceptions in the non-bonded and\n  ghost-14 forces\n" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::toString
        
            typedef ::QString ( ::SireOpenMM::PerturbableOpenMMMolecule::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::toString );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::typeName );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireOpenMM::PerturbableOpenMMMolecule::what
        
            typedef char const * ( ::SireOpenMM::PerturbableOpenMMMolecule::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireOpenMM::PerturbableOpenMMMolecule::what );
            
            PerturbableOpenMMMolecule_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        PerturbableOpenMMMolecule_exposer.staticmethod( "typeName" );
        PerturbableOpenMMMolecule_exposer.def( "__copy__", &__copy__);
        PerturbableOpenMMMolecule_exposer.def( "__deepcopy__", &__copy__);
        PerturbableOpenMMMolecule_exposer.def( "clone", &__copy__);
        PerturbableOpenMMMolecule_exposer.def( "__str__", &__str__< ::SireOpenMM::PerturbableOpenMMMolecule > );
        PerturbableOpenMMMolecule_exposer.def( "__repr__", &__str__< ::SireOpenMM::PerturbableOpenMMMolecule > );
        PerturbableOpenMMMolecule_exposer.def( "__str__", &__str__< ::SireOpenMM::PerturbableOpenMMMolecule > );
        PerturbableOpenMMMolecule_exposer.def( "__repr__", &__str__< ::SireOpenMM::PerturbableOpenMMMolecule > );
    }

}
