// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "SupraSubSimPacket.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "suprasubsimpacket.h"

#include <QDebug>

#include "suprasubsimpacket.h"

SireMove::SupraSubSimPacket __copy__(const SireMove::SupraSubSimPacket &other){ return SireMove::SupraSubSimPacket(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMove::SupraSubSimPacket&){ return "SireMove::SupraSubSimPacket";}

void register_SupraSubSimPacket_class(){

    { //::SireMove::SupraSubSimPacket
        typedef bp::class_< SireMove::SupraSubSimPacket, bp::bases< SireCluster::WorkPacketBase > > SupraSubSimPacket_exposer_t;
        SupraSubSimPacket_exposer_t SupraSubSimPacket_exposer = SupraSubSimPacket_exposer_t( "SupraSubSimPacket", "This is a workpacket that is used to run part of a SupraSubSim simulation\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope SupraSubSimPacket_scope( SupraSubSimPacket_exposer );
        SupraSubSimPacket_exposer.def( bp::init< SireMove::SupraSubSystem const &, SireMove::SupraSubMoves const &, int, bp::optional< bool > >(( bp::arg("system"), bp::arg("moves"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) ), "Construct a work packet to perform nmoves sub-moves (in moves) on\nthe sub-system system, recording statistics if record_stats is true") );
        SupraSubSimPacket_exposer.def( bp::init< SireMove::SupraSubSimPacket const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::SupraSubSimPacket::approximatePacketSize
        
            typedef int ( ::SireMove::SupraSubSimPacket::*approximatePacketSize_function_type)(  ) const;
            approximatePacketSize_function_type approximatePacketSize_function_value( &::SireMove::SupraSubSimPacket::approximatePacketSize );
            
            SupraSubSimPacket_exposer.def( 
                "approximatePacketSize"
                , approximatePacketSize_function_value
                , "This is likely to be large" );
        
        }
        { //::SireMove::SupraSubSimPacket::hasFinished
        
            typedef bool ( ::SireMove::SupraSubSimPacket::*hasFinished_function_type)(  ) const;
            hasFinished_function_type hasFinished_function_value( &::SireMove::SupraSubSimPacket::hasFinished );
            
            SupraSubSimPacket_exposer.def( 
                "hasFinished"
                , hasFinished_function_value
                , "Return whether or not this work packet has finished" );
        
        }
        { //::SireMove::SupraSubSimPacket::nSubCompleted
        
            typedef int ( ::SireMove::SupraSubSimPacket::*nSubCompleted_function_type)(  ) const;
            nSubCompleted_function_type nSubCompleted_function_value( &::SireMove::SupraSubSimPacket::nSubCompleted );
            
            SupraSubSimPacket_exposer.def( 
                "nSubCompleted"
                , nSubCompleted_function_value
                , "Return the number of completed sub-moves" );
        
        }
        { //::SireMove::SupraSubSimPacket::nSubMoves
        
            typedef int ( ::SireMove::SupraSubSimPacket::*nSubMoves_function_type)(  ) const;
            nSubMoves_function_type nSubMoves_function_value( &::SireMove::SupraSubSimPacket::nSubMoves );
            
            SupraSubSimPacket_exposer.def( 
                "nSubMoves"
                , nSubMoves_function_value
                , "Return the number of sub-moves to be applied to the sub-system" );
        
        }
        SupraSubSimPacket_exposer.def( bp::self != bp::self );
        { //::SireMove::SupraSubSimPacket::operator=
        
            typedef ::SireMove::SupraSubSimPacket & ( ::SireMove::SupraSubSimPacket::*assign_function_type)( ::SireMove::SupraSubSimPacket const & ) ;
            assign_function_type assign_function_value( &::SireMove::SupraSubSimPacket::operator= );
            
            SupraSubSimPacket_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        SupraSubSimPacket_exposer.def( bp::self == bp::self );
        { //::SireMove::SupraSubSimPacket::recordingSubStatistics
        
            typedef bool ( ::SireMove::SupraSubSimPacket::*recordingSubStatistics_function_type)(  ) const;
            recordingSubStatistics_function_type recordingSubStatistics_function_value( &::SireMove::SupraSubSimPacket::recordingSubStatistics );
            
            SupraSubSimPacket_exposer.def( 
                "recordingSubStatistics"
                , recordingSubStatistics_function_value
                , "Return whether or not we are recording statistics during the sub-moves" );
        
        }
        { //::SireMove::SupraSubSimPacket::shouldPack
        
            typedef bool ( ::SireMove::SupraSubSimPacket::*shouldPack_function_type)(  ) const;
            shouldPack_function_type shouldPack_function_value( &::SireMove::SupraSubSimPacket::shouldPack );
            
            SupraSubSimPacket_exposer.def( 
                "shouldPack"
                , shouldPack_function_value
                , "We probably shouldnt pack this workpacket as it is already\nheavily packed, and if it is large, then it is likely already\npacked to disk" );
        
        }
        { //::SireMove::SupraSubSimPacket::subMoves
        
            typedef ::SireMove::SupraSubMoves const & ( ::SireMove::SupraSubSimPacket::*subMoves_function_type)(  ) const;
            subMoves_function_type subMoves_function_value( &::SireMove::SupraSubSimPacket::subMoves );
            
            SupraSubSimPacket_exposer.def( 
                "subMoves"
                , subMoves_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the moves" );
        
        }
        { //::SireMove::SupraSubSimPacket::subSystem
        
            typedef ::SireMove::SupraSubSystem const & ( ::SireMove::SupraSubSimPacket::*subSystem_function_type)(  ) const;
            subSystem_function_type subSystem_function_value( &::SireMove::SupraSubSimPacket::subSystem );
            
            SupraSubSimPacket_exposer.def( 
                "subSystem"
                , subSystem_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the sub-system" );
        
        }
        { //::SireMove::SupraSubSimPacket::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::SupraSubSimPacket::typeName );
            
            SupraSubSimPacket_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMove::SupraSubSimPacket::what
        
            typedef char const * ( ::SireMove::SupraSubSimPacket::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMove::SupraSubSimPacket::what );
            
            SupraSubSimPacket_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        SupraSubSimPacket_exposer.staticmethod( "typeName" );
        SupraSubSimPacket_exposer.def( "__copy__", &__copy__);
        SupraSubSimPacket_exposer.def( "__deepcopy__", &__copy__);
        SupraSubSimPacket_exposer.def( "clone", &__copy__);
        SupraSubSimPacket_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::SupraSubSimPacket >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SupraSubSimPacket_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::SupraSubSimPacket >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SupraSubSimPacket_exposer.def( "__setstate__", &__setstate__base64< ::SireMove::SupraSubSimPacket > );
        SupraSubSimPacket_exposer.def( "__getstate__", &__getstate__base64< ::SireMove::SupraSubSimPacket > );
        SupraSubSimPacket_exposer.def( "__str__", &pvt_get_name);
        SupraSubSimPacket_exposer.def( "__repr__", &pvt_get_name);
    }

}
