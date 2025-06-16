library(data.table)
library(optparse)
library(Rlabkey)
library(dplyr)


option_list <- list(
  make_option(c("-d", "--disease"), type = "character", default = NULL,
              help = "Disease term by which to define a row as case or control by (e.g., 'hypertension')- see readme for how to specify", metavar = "character"),
  make_option(c("-l", "--list"), type = "character", default = NULL,
              help = "file with list of patients (in the first column) that need to be annotated- output annotations will be merged to this file (including other columns)- if none is supplied all patients in 'panels_applied' will be annotated", metavar = "character"),
  make_option(c("-m", "--merge_by"), type = "character", default = "plate_key",
              help = "If the supplied --list input file needs to be searched according to a column other than the first one- please input that column name here", 
              metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "out.tsv",
              help = "output file into which the data will be put- default is 'out.tsv'", metavar = "character"),
  make_option(c("-c", "--control_age_cutoff"), type = "character", default = "out.tsv",
              help = "Year of birth at which to remove controls- defined according to '--disease'", metavar = "character"),
  make_option(c("-p", "--get_paths"), type = "character", default = "TRUE",
              help = "option to get paths to cram files and genome versions- default = TRUE", metavar = "character"),
  make_option(c("-n", "--number"), type = "integer", default = NULL,
              help = "option to take first n rows in from '--list' option", metavar = "character"),
  make_option(c("-a", "--annotate"), type = "character", default = "TRUE",
              help = "option to annotate rows with disease according to 'panels applied, phenotype and hpo terms'", metavar = "character"),
  make_option(c("-e", "--cancer_genomes"), type = "character", default = NULL,
              help = "include cancer genomes as controls", metavar = "character"),
  make_option(c("-k", "--order_by"), type = "character", default = NULL,
              help = "order the input and output file according to a particular column- specify column name", metavar = "character"),
  make_option(c("-s", "--suffix"), type = "character", default = NULL,
              help = "remove suffix of column on which search is being performed", metavar = "character"),
  make_option(c("-r", "--prefix"), type = "character", default = NULL,
              help = "remove prefix of column on which search is being performed", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



query=function(patients){
     # Create SQL query that joins the required tables
    phenotype_columns_100k=paste0("pa.rare_diseases_family_id,  pa.panel_name,  pa.phenotype,  rd.normalised_hpo_term,  rd.hpo_term")
    sql_columns_100k <- paste0("pa.sample_id AS plate_key, pa.participant_id, p.participant_phenotyped_sex AS gender,  p.yob AS year_of_birth,  p.genetically_inferred_ancestry_thr AS ancestry")
    tables_100k=paste0(" FROM panels_applied pa LEFT JOIN rare_diseases_participant_phenotype rd   ON pa.participant_id = rd.participant_id LEFT JOIN participant_summary p   ON pa.participant_id = p.participant_id")
    phenotype_columns_gms=paste0("pa.referral_id,  pa.panel_name,  pa.phenotype,  obs.normalised_hpo_term,  obs.code_description")
    sql_columns_gms <- paste0(" pa.platekey,  pa.participant_id, p.administrative_gender AS gender,  p.participant_year_of_birth AS year_of_birth,  p.ethnicity_description AS ancestry")
    tables_gms=paste0(" FROM panels_applied pa LEFT JOIN observation obs   ON pa.participant_id = obs.participant_id LEFT JOIN participant p   ON pa.participant_id = p.participant_id")

    sql_str_100k=paste0(sql_columns_100k,tables_100k)
    sql_str_gms=paste0(sql_columns_gms,tables_gms)
    if(opt$annotate=="TRUE"){
        sql_str_100k=paste0(phenotype_columns_100k,', ',sql_str_100k)
        sql_str_gms=paste0(phenotype_columns_gms,', ',sql_str_gms)
    }

    if(opt$get_paths=='TRUE'){
        sql_100k_paths=paste0('path.platekey,path.genome_build,path.file_path,')
        sql_str_100k=paste(sql_100k_paths,sql_str_100k,'LEFT JOIN genome_file_paths_and_types path ON pa.participant_id = path.participant_id')
        sql_gms_paths=paste0('path.platekey,path.genome_build,path.path,')
        sql_str_gms=paste(sql_gms_paths,sql_str_gms,'LEFT JOIN genome_file_paths_and_types path ON pa.participant_id = path.participant_id')
        sql_str_100k=paste0(sql_str_100k," WHERE path.file_sub_type = 'BAM'")
        sql_str_gms=paste0(sql_str_gms," WHERE path.file_sub_type = 'CRAM'")
    }

    search_gel=TRUE
    search_gms=TRUE
    bind=FALSE
    # if(!is.null(opt$list)){
        # if(file.exists(opt$list)){
            if(opt$merge_by=='participant_id'){search_pattern='pp'}else{search_pattern='LP5'}
            n_platekey_gms=sum(grepl(pattern=search_pattern,patients),na.rm=T)/length(patients)
            print(n_platekey_gms)
            if((n_platekey_gms<1)){
                search_gel=TRUE
            }else{
                search_gel=FALSE
            }
            if((n_platekey_gms>0)){
                search_gms=TRUE
            }else{
                search_gms=FALSE
            }
            
            if(opt$merge_by=='plate_key'){gel_query_col='sample_id';gms_query_col='platekey'}else{gel_query_col=opt$merge_by;gms_query_col=opt$merge_by}
            selection1=paste0("pa.",gel_query_col," IN ('", paste(patients,collapse="\',\'"), "');")
            selection2=paste0("pa.",gms_query_col," IN ('", paste(patients,collapse="\',\'"), "');")
            if(grepl(pattern='WHERE',sql_str_100k)){
                selection1=paste0("AND ",selection1)
                selection2=paste0("AND ",selection2)
            }else{
                selection1=paste0("WHERE ",selection1)
                selection2=paste0("WHERE ",selection2)
            }
            bind=TRUE
        # }else{
        #     if(opt$merge_by=='plate_key'){gel_query_col='sample_id'}else{gel_query_col=opt$merge_by}
        #     selection1=paste0("WHERE pa.",gel_query_col," IN ('", opt$list, "');")
        #     selection2=paste0("WHERE pa.platekey IN ('", opt$list, "');")
        # }
        sql_str_100k=paste(sql_str_100k,selection1)
        sql_str_gms=paste(sql_str_gms,selection2)
    # }else{
    #     sql_str_100k=paste0(sql_str_100k,';')
    #     sql_str_gms=paste0(sql_str_gms,';')
    # }


    sql_str_100k=paste("SELECT",sql_str_100k)
    sql_str_gms=paste("SELECT",sql_str_gms)
    # print(sql_str_100k)
    # print(sql_str_gms)


    if(search_gel){
        gel_data <- unique(labkey.executeSql(
            baseUrl="https://labkey-embassy.gel.zone/labkey", 
            folderPath="/main-programme/main-programme_v19_2024-10-31", 
            schemaName="lists", 
            containerFilter=NULL, 
            sql = sql_str_100k
        ))
        if(sum(c("Plate Key","Platekey")%in%colnames(gel_data))==2){gel_data=gel_data[,colnames(gel_data!='Plake Key')]}
    }else{gel_data=data.frame()}
    # print(gel_data)
    # print(head(gel_data))
    # Execute the SQL
    # print('1')
    if(search_gms){
        gms_data <- unique(labkey.executeSql(
            baseUrl = "https://labkey-embassy.gel.zone/labkey",
            folderPath = "/nhs-gms/nhs-gms-release_v4_2024-08-22",
            schemaName = "lists",
            containerFilter = NULL,
            sql = sql_str_gms
        ))
        if(sum(c("Plate Key","Platekey")%in%colnames(gms_data))==2){gms_data=gms_data[,colnames(gms_data!='Plake Key')]}
    }else{gms_data=data.frame()}
    
    # print('2')
    if(search_gel | search_gms){
    column_translations=data.frame(real_names=c("Plate Key","Platekey",'File Path','file_path',"Genome Build","Participant Id","Rare Diseases Family Id","Referral Id","Panel Name","Phenotype","Normalised Hpo Term","Hpo Term","Code Description","Gender","Year Of Birth","Ancestry"),
                                translations=c("plate_key","plate_key",'path','path',"genome_build","participant_id","family_id","family_id","panel_name","phenotype","normalised_hpo_term","code_description","code_description","gender","year_of_birth","ancestry"),
                        stringsAsFactors=F)

    for(col in column_translations$real_names){
        translation=column_translations$translations[column_translations$real_names==col]
        if(search_gel){
            colnames(gel_data)[colnames(gel_data)==col]=translation
        }
        if(search_gms){
            colnames(gms_data)[colnames(gms_data)==col]=translation
        }
    }

    if(search_gel & search_gms){
        if(nrow(gel_data)>0 | nrow(gms_data)>0){
            data_out=rbind(gel_data,gms_data)
        }
    }else{
        if(search_gel){
            if(nrow(gel_data)>0){
                data_out=gel_data
            }
        }
        if(search_gms){
            if(nrow(gms_data)>0){
                data_out=gms_data
            }
        }
    }
# print('3')
    if(sum(colnames(data_out)=='plate_key')>1){idx=which(colnames(data_out)=='plate_key')[2];data_out=data_out[,-idx]}

    if(!is.null(opt$disease)){
        data_out$case_or_control='control'
        # data_out$case_or_control[grep(pattern=opt$disease,data_out$collapsed_panel_name)]='case'
        # data_out$case_or_control[is.na(data_out$case_or_control)]='control'
        filters=unlist(strsplit(opt$disease,','))
        for(filter in filters){
            filter=tolower(filter)
            if(grepl(pattern='=|equal', filter)){
                params=unlist(strsplit(filter,"=|equal"))
                column_oi=gsub(" ","",params[1])
                value=gsub("^\\W+","",params[2])
                data_out[data_out[,column_oi]==value,'case_or_control']='case'
            }else if(grepl(pattern='contains', filter)){
                params=unlist(strsplit(filter,"contains"))
                column_oi=gsub(" ","",params[1])
                value=gsub("^\\W+","",params[2])
                data_out[grepl(value,data_out[,column_oi]),'case_or_control']='case'
            }
        }
        if(!is.null(opt$control_age_cuttoff)){
            data_out=data_out[!(data_out$year_of_birth>opt$control_age_cuttoff & data_out$case_or_control=='control'),]
        }
    }
    # print('4')
    # print(class(data_out))
    if(opt$annotate=="TRUE"){
        if("genome_build"%in%colnames(data_out)){
            data_out=data_out%>%group_by(plate_key,participant_id,path,genome_build)%>%
                mutate(gender=unique(gender),year_of_birth=unique(year_of_birth),ancestry=unique(ancestry),
                    panel_name_collapsed=paste(unique(panel_name),collapse=';'),collapsed_phenotype=paste(unique(phenotype),collapse=';'),
                    collapsed_hpo_term=paste(unique(normalised_hpo_term),collapse=';'),
                    collapsed_hpo_description=paste(unique(code_description),collapse=';'))%>%
                    select(plate_key,participant_id,path,genome_build,participant_id,gender,year_of_birth,ancestry,panel_name_collapsed,collapsed_phenotype,collapsed_hpo_term,collapsed_hpo_description,case_or_control)%>%
                unique()
        }else{
        data_out=data_out%>%group_by(plate_key,participant_id)%>%
            mutate(gender=unique(gender),year_of_birth=unique(year_of_birth),ancestry=unique(ancestry),
                panel_name_collapsed=paste(unique(panel_name),collapse=';'),collapsed_phenotype=paste(unique(phenotype),collapse=';'),
                collapsed_hpo_term=paste(unique(normalised_hpo_term),collapse=';'),
                collapsed_hpo_description=paste(unique(code_description),collapse=';'))%>%
                select(plate_key,participant_id,gender,year_of_birth,ancestry,panel_name_collapsed,collapsed_phenotype,collapsed_hpo_term,collapsed_hpo_description,case_or_control)%>%
            unique()
        }
    }
    # print('5')
    # print(data_out)
    # print(dim(data_out))
    data_out=data.frame(data_out,stringsAsFactors=F)
    
    if(bind){
        data_out=merge(data,data_out,by=opt$merge_by)
    }
    return(data_out)
}
}


if(!is.null(opt$list)){
    data=fread(opt$list,sep='\t',fill=T)
    data=data.frame(data,stringsAsFactors=F)
    if(!is.null(opt$order_by)){
        data=data[order(data[,opt$order_by],decreasing=T),]
        # print(data)
    }
    if(!is.null(opt$number)){
        data=head(data,n=opt$n)
    }
    colnames(data)[1]=opt$merge_by
    patients=data.frame(data,stringsAsFactors=F)[,opt$merge_by]
    bind=TRUE
}else{
    patients_gel=labkey.selectRows(
        baseUrl="https://labkey-embassy.gel.zone/labkey", 
        folderPath="/main-programme/main-programme_v19_2024-10-31", 
        schemaName="lists", 
        queryName="panels_applied", 
        viewName="", 
        colSelect="sample_id", 
        colFilter=NULL, 
        containerFilter=NULL, 
        colNameOpt="rname"
    )
    colnames(patients_gel)[1]='plate_key'
    patients_gms=labkey.selectRows(
        baseUrl="https://labkey-embassy.gel.zone/labkey", 
        folderPath="/nhs-gms/nhs-gms-release_v4_2024-08-22", 
        schemaName="lists", 
        queryName="panels_applied", 
        viewName="", 
        colSelect="platekey", 
        colFilter=NULL, 
        containerFilter=NULL, 
        colNameOpt="rname"
    )
    colnames(patients_gms)[1]='plate_key'
    data=rbind(patients_gel,patients_gms)
    patients=data$plate_key
    if(!is.null(opt$number)){
        patients=sample(patients,opt$number)
        print(patients)
    }
}
print(length(patients))
if(length(patients)<1000){
    data_out=query(patients)
}else{
    patient_chunks=split(patients,ceiling(seq_along(patients)/1000))
    all_data=lapply(patient_chunks,function(patients){
        return(query(patients))
    })
    data_out=do.call('rbind',all_data)
}

# print('4')
# print(head(data_out))

if(!is.null(opt$cancer_genomes)){
    gel_cancer_genomes <- labkey.selectRows(
        baseUrl="https://labkey-embassy.gel.zone/labkey", 
        folderPath="/main-programme/main-programme_v19_2024-10-31", 
        schemaName="lists", 
        queryName="cancer_analysis", 
        viewName="",
        colSelect="participant_id,germline_sample_platekey,sex,year_of_birth,germline_bam", 
        colFilter=NULL, 
        containerFilter=NULL, 
        colNameOpt="rname"
    )
    colnames(gel_cancer_genomes)=c('participant_id','plate_key','gender','year_of_birth','path')
    gel_cancer_genomes$genome_build='GRCh38'
    cols_not_present=colnames(data_out)[colnames(data_out)%in%colnames(gel_cancer_genomes)]
    for(col in cols_not_present){
        gel_cancer_genomes[,col]=NA
    }
    gel_cancer_genomes=gel_cancer_genomes[,colnames(data_out)]

    gms_cancer_genomes=labkey.selectRows(
        baseUrl="https://labkey-embassy.gel.zone/labkey", 
        folderPath="/nhs-gms/nhs-gms-release_v4_2024-08-22", 
        schemaName="lists", 
        queryName="cancer_analysis", 
        viewName="", 
        colSelect="participant_id,germline_sample_platekey,sex,year_of_birth,germline_alignment_path", 
        colFilter=NULL, 
        containerFilter=NULL, 
        colNameOpt="rname"
    )
    gms_cancer_genomes$genome_build='GRCh38'
    cols_not_present=colnames(data_out)[colnames(data_out)%in%colnames(gms_cancer_genomes)]
    for(col in cols_not_present){
        gms_cancer_genomes[,col]=NA
    }
    gms_cancer_genomes=gms_cancer_genomes[,colnames(data_out)]
    data_out=rbind(data_out,gel_cancer_genomes,gms_cancer_genomes)
}

datagr37=data_out[data_out$genome_build=='GRCh37',]
datagr37=datagr37[!(datagr37$plate_key%in%data_out$plate_key),]
if(nrow(datagr37)>0){
    data_out=rbind(data_out,datagr37)
}
if(!is.null(opt$order_by)){
    data_out=data_out[order(data_out[,opt$order_by],decreasing=T),]
}
write.table(data_out,opt$output,sep='\t',row.names=F,quote=F)

# }else{
#     print('no recognisable platekeys/ participant ids in query')
#     print(patients)
#     print(search_gel)
#     print(search_gms)
# }