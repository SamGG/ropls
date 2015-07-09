\name{foods}

\alias{foods}

\docType{data}

\title{
  Food consumption patterns accross European countries (FOODS)
}

\description{
  The relative consumption of 20 food items was compiled for 16
  countries. The values range between 0 and 100 percent and a high value
  corresponds to a high consumption. The dataset contains 3 missing data.
}

\usage{
data(foods)
}

\format{
  A data frame with the following parameters:
  \describe{
    \item{\code{Country}}{
      Name of the country
    }
    \item{\code{Gr_Coffe}}{
      Ground Coffee
    }
    \item{\code{Inst_Coffe}}{
      Instant Coffee
    }
    \item{\code{Tea}}{
      Tea
    }
    \item{\code{Sweetner}}{
      Sweetner
    }
    \item{\code{Biscuits}}{
      Biscuits
    }
    \item{\code{Pa_Soup}}{
      Powder Soup
    }
    \item{\code{Ti_Soup}}{
      Tin Soup
    }
    \item{\code{In_Potat}}{
      Instant Potatoes
    }
    \item{\code{Fro_Fish}}{
      Frozen Fish
    }
    \item{\code{Fro_Veg}}{
      Frozen Vegetables
    }
    \item{\code{Apples}}{
     Apples
    }
    \item{\code{Oranges}}{
      Oranges
    }
    \item{\code{Ti_Fruit}}{
      Tin Fruit
    }
    \item{\code{Jam}}{
      Jam
    }
    \item{\code{Garlic}}{
      Garlic
    }
    \item{\code{Butter}}{
      Butter
    }
    \item{\code{Margarine}}{
      Margarine
    }
    \item{\code{Olive_Oil}}{
      Olive Oil
    }
    \item{\code{Yoghurt}}{
      Yoghurt
    }
    \item{\code{Crisp_Brea}}{
      Crisp Bread
    }
  }
}

\value{
  Data frame (numeric type except the first column, which can be
  transformed into row names) with 16 rows and 21 columns, corresponding to the 'Country' (column 1), followed by the consumption
  of each of the 20 food items (columns 2-21). For details see the 'Format'
  section above.
}

\references{
Eriksson et al. (2006). Multi- and Megarvariate Data Analysis. Umetrics
Academy. pp.10, 33, 48.
}

\keyword{datasets}
