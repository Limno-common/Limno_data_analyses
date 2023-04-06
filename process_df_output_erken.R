# Function to write the DF raw output file into csv tables

# Update 2023-04-06: addition of row_id

library(data.table)
library(xml2)
library(stringr)

process_df_output_erken = function(df_txt_output, folder_out = "."){
  # There are three types of measurements: composition, biomass, and pi-curve
  # Write a table for each of these
  df_composition = data.table()
  df_biomass = data.table()
  df_picurve = data.table()
  
  # Read file
  lines = readLines(df_txt_output)
  
  # Locate the start indices and end indices of the readings
  start_pattern = "<measurement type="
  start_ind = str_which(lines, start_pattern)
  end_ind = c(start_ind[-1L] - 1L, length(lines))
  
  # Read entries
  row_id_comp = 1L
  row_id_bio = 1L
  row_id_pic = 1L
  for(i in seq_along(start_ind)){
    if(grepl("type=\"composition\"", lines[start_ind[i]])){
      the_df = "df_composition"
      the_row_id = "row_id_comp"
    }else if(grepl("type=\"biomass\"", lines[start_ind[i]])){
      the_df = "df_biomass"
      the_row_id = "row_id_bio"
    }else if(grepl("type=\"pi_curve\"", lines[start_ind[i]])){
      the_df = "df_picurve"
      the_row_id = "row_id_pic"
    }else{
      stop("Type in line ", start_ind[i], " is not one of composition, biomass, or pi_curve")
    }
    
    # Extract the xml data entry
    measurement_string = paste(lines[start_ind[i]:end_ind[i]], collapse = "\n")
    m_start_ind = str_locate(measurement_string, "<measurement")
    m_end_ind = str_locate(measurement_string, "</measurement>")
    measurement = substring(measurement_string, first = m_start_ind[1L], last = m_end_ind[2L])
    
    xml_file = read_xml(measurement)
    xml = as_list(xml_file)
    
    # Extract the date, darkrate, and darkratevar
    the_date = attr(xml[["measurement"]], "date")
    the_darkrate = attr(xml[["measurement"]], "darkrate")
    the_darkratevar = attr(xml[["measurement"]], "darkratevar")
    
    # Get parameter names and values
    lst_to_add = list(row_id = get(the_row_id),
                      date = the_date,
                      darkrate = the_darkrate,
                      darkratevar = the_darkratevar)
    if(the_df == "df_composition"){
      for(j in seq_len(length(xml[["measurement"]]))){
        col_name = paste(attr(xml[["measurement"]][j][[1L]], "name"),
                         attr(xml[["measurement"]][j][[1L]], "wavelength"),
                         sep = "_")
        the_vals = xml[["measurement"]][j][[1L]][[1L]]
        if(the_vals == " ") next
        lst_to_add[[col_name]] = as.numeric(str_split(the_vals, " ")[[1L]])
      }
    }else if(the_df == "df_biomass"){
      the_vals = xml[["measurement"]][["throughflow"]][[1L]]
      lst_to_add[["throughflow"]] = as.numeric(str_split(the_vals, " ")[[1L]])
      the_vals = xml[["measurement"]][["decay"]][[1L]]
      lst_to_add[["decay"]] = as.numeric(str_split(the_vals, " ")[[1L]])
      
      # Throughflow count is shorter than decay, so add NAs
      lst_to_add[["throughflow"]] = c(lst_to_add[["throughflow"]],
                                      rep(NA, length(lst_to_add[["decay"]]) -
                                            length(lst_to_add[["throughflow"]])))
      
    }else if(the_df == "df_picurve"){
      for(j in seq_len(length(xml[["measurement"]]))){
        col_name = paste0("light_", attr(xml[["measurement"]][j][[1L]], "light"))
        the_vals = xml[["measurement"]][j][[1L]][[1L]]
        if(the_vals == " ") next
        lst_to_add[[col_name]] = as.numeric(str_split(the_vals, " ")[[1L]])
      }
    }
    
    assign(the_df, rbindlist(list(get(the_df), lst_to_add), fill = T))
    assign(the_row_id, get(the_row_id) + 1L)
  }
  
  # Sort the tables by date
  df_composition[, date := as.POSIXct(date, tz = "UTC")]
  df_biomass[, date := as.POSIXct(date, tz = "UTC")]
  df_picurve[, date := as.POSIXct(date, tz = "UTC")]
  setorder(df_composition, date)
  setorder(df_biomass, date)
  setorder(df_picurve, date)
  
  # Write tables
  fwrite(df_composition, file.path(folder_out, "DF_output_composition.csv"), dateTimeAs = "write.csv")
  fwrite(df_biomass, file.path(folder_out, "DF_output_biomass.csv"), dateTimeAs = "write.csv")
  fwrite(df_picurve, file.path(folder_out, "DF_output_picurve.csv"), dateTimeAs = "write.csv")
}
