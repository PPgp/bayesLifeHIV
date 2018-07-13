
hiv.countries.pred <- function(meta)
    return(data.frame(country_code=meta$regions$country_code[meta$regions$hiv.pred],
                 country_name=meta$regions$country_name[meta$regions$hiv.pred])
            )

hiv.countries.est <- function(meta)
    return(data.frame(country_code=meta$regions$country_code[meta$regions$hiv.est],
                      country_name=meta$regions$country_name[meta$regions$hiv.est])
    )