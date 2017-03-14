OTB_CREATE_APPLICATION(NAME MultiscaleSegmentation
                       SOURCES MultiscaleSegmentation.cxx otbWaveletGenerator.cxx
                       LINK_LIBRARIES OTBIO OTBCommon OTBBasicFilters OTBOGRAdapters)
