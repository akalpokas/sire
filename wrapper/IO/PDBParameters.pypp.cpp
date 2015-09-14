// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PDBParameters.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/stringmangler.h"

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomeditor.h"

#include "SireMol/atomelements.h"

#include "SireMol/cgeditor.h"

#include "SireMol/chain.h"

#include "SireMol/chaineditor.h"

#include "SireMol/connectivity.h"

#include "SireMol/cuttingfunction.h"

#include "SireMol/element.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleditor.h"

#include "SireMol/molidentifier.h"

#include "SireMol/molidx.h"

#include "SireMol/mover.hpp"

#include "SireMol/reseditor.h"

#include "SireMol/residue.h"

#include "SireMol/segeditor.h"

#include "SireMol/segment.h"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireStream/sharestrings.h"

#include "SireUnits/units.h"

#include "pdb.h"

#include "tostring.h"

#include <QByteArray>

#include <QDebug>

#include <QMap>

#include <QString>

#include <QStringList>

#include <QTextStream>

#include "pdb.h"

SireIO::PDBParameters __copy__(const SireIO::PDBParameters &other){ return SireIO::PDBParameters(other); }

const char* pvt_get_name(const SireIO::PDBParameters&){ return "SireIO::PDBParameters";}

void register_PDBParameters_class(){

    { //::SireIO::PDBParameters
        typedef bp::class_< SireIO::PDBParameters, bp::bases< SireIO::IOParametersBase > > PDBParameters_exposer_t;
        PDBParameters_exposer_t PDBParameters_exposer = PDBParameters_exposer_t( "PDBParameters", bp::init< >() );
        bp::scope PDBParameters_scope( PDBParameters_exposer );
        { //::SireIO::PDBParameters::alternatives
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*alternatives_function_type)(  ) const;
            alternatives_function_type alternatives_function_value( &::SireIO::PDBParameters::alternatives );
            
            PDBParameters_exposer.def( 
                "alternatives"
                , alternatives_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::animationFrameSelector
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*animationFrameSelector_function_type)(  ) const;
            animationFrameSelector_function_type animationFrameSelector_function_value( &::SireIO::PDBParameters::animationFrameSelector );
            
            PDBParameters_exposer.def( 
                "animationFrameSelector"
                , animationFrameSelector_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::animationFrames
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*animationFrames_function_type)(  ) const;
            animationFrames_function_type animationFrames_function_value( &::SireIO::PDBParameters::animationFrames );
            
            PDBParameters_exposer.def( 
                "animationFrames"
                , animationFrames_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::atomNameMangler
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*atomNameMangler_function_type)(  ) const;
            atomNameMangler_function_type atomNameMangler_function_value( &::SireIO::PDBParameters::atomNameMangler );
            
            PDBParameters_exposer.def( 
                "atomNameMangler"
                , atomNameMangler_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::bFactor
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*bFactor_function_type)(  ) const;
            bFactor_function_type bFactor_function_value( &::SireIO::PDBParameters::bFactor );
            
            PDBParameters_exposer.def( 
                "bFactor"
                , bFactor_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::chainNameMangler
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*chainNameMangler_function_type)(  ) const;
            chainNameMangler_function_type chainNameMangler_function_value( &::SireIO::PDBParameters::chainNameMangler );
            
            PDBParameters_exposer.def( 
                "chainNameMangler"
                , chainNameMangler_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::formalCharge
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*formalCharge_function_type)(  ) const;
            formalCharge_function_type formalCharge_function_value( &::SireIO::PDBParameters::formalCharge );
            
            PDBParameters_exposer.def( 
                "formalCharge"
                , formalCharge_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::iCode
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*iCode_function_type)(  ) const;
            iCode_function_type iCode_function_value( &::SireIO::PDBParameters::iCode );
            
            PDBParameters_exposer.def( 
                "iCode"
                , iCode_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::pdbAtomName
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*pdbAtomName_function_type)(  ) const;
            pdbAtomName_function_type pdbAtomName_function_value( &::SireIO::PDBParameters::pdbAtomName );
            
            PDBParameters_exposer.def( 
                "pdbAtomName"
                , pdbAtomName_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::pdbChainName
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*pdbChainName_function_type)(  ) const;
            pdbChainName_function_type pdbChainName_function_value( &::SireIO::PDBParameters::pdbChainName );
            
            PDBParameters_exposer.def( 
                "pdbChainName"
                , pdbChainName_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::pdbResidueName
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*pdbResidueName_function_type)(  ) const;
            pdbResidueName_function_type pdbResidueName_function_value( &::SireIO::PDBParameters::pdbResidueName );
            
            PDBParameters_exposer.def( 
                "pdbResidueName"
                , pdbResidueName_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::pdbSegmentName
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*pdbSegmentName_function_type)(  ) const;
            pdbSegmentName_function_type pdbSegmentName_function_value( &::SireIO::PDBParameters::pdbSegmentName );
            
            PDBParameters_exposer.def( 
                "pdbSegmentName"
                , pdbSegmentName_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::residueNameMangler
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*residueNameMangler_function_type)(  ) const;
            residueNameMangler_function_type residueNameMangler_function_value( &::SireIO::PDBParameters::residueNameMangler );
            
            PDBParameters_exposer.def( 
                "residueNameMangler"
                , residueNameMangler_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireIO::PDBParameters::segmentNameMangler
        
            typedef ::SireBase::PropertyName const & ( ::SireIO::PDBParameters::*segmentNameMangler_function_type)(  ) const;
            segmentNameMangler_function_type segmentNameMangler_function_value( &::SireIO::PDBParameters::segmentNameMangler );
            
            PDBParameters_exposer.def( 
                "segmentNameMangler"
                , segmentNameMangler_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        PDBParameters_exposer.def( "__copy__", &__copy__);
        PDBParameters_exposer.def( "__deepcopy__", &__copy__);
        PDBParameters_exposer.def( "clone", &__copy__);
        PDBParameters_exposer.def( "__str__", &pvt_get_name);
        PDBParameters_exposer.def( "__repr__", &pvt_get_name);
    }

}
