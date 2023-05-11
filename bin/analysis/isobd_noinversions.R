#relatedness analysis
library(tidyverse)
library(data.table)
library(geosphere)
library(adegenet)

mrel = fread('~/Projects/td_je_angam_2022/data/relatedness/m_form.res')
srel = fread('~/Projects/td_je_angam_2022/data/relatedness/s_form.res')
mdata = fread('~/Projects/td_je_angam_2022/metadata/metadata_with_insecticide_info.csv')



m_data <- mdata[mdata$Form == 'M']
m_lalo <- select(m_data, Lat, Long)
m_lalo$id <- seq(0, 158)
aj = left_join(mrel, m_lalo, by=c('a' = 'id'))
bj = left_join(aj, m_lalo, by=c('b' = 'id'))
bj$dist = distVincentyEllipsoid(bj[,c('Lat.x','Long.x')], bj[,c('Lat.y','Long.y')])
x = bj %>% select(a,b, rab) %>% pivot_wider(names_from = a, values_from = rab)
x = x[,-1]
y = bj %>% select(a,b, dist) %>% pivot_wider(names_from = a, values_from = dist)
y = y[,-1]

s_data <- mdata[mdata$Form == 'S']
s_lalo <- select(s_data, Lat, Long)
s_lalo$id <- seq(0, 154)
sj = left_join(srel, s_lalo, by=c('a' = 'id'))
cj = left_join(sj, s_lalo, by=c('b' = 'id'))
cj$dist = distVincentyEllipsoid(cj[,c('Lat.x','Long.x')], cj[,c('Lat.y','Long.y')])
t = cj %>% select(a,b, rab) %>% pivot_wider(names_from = a, values_from = rab)
t = x[,-1]
z = cj %>% select(a,b, dist) %>% pivot_wider(names_from = a, values_from = dist)
z = y[,-1]


misobd = ggplot(bj, aes(x=dist,y=KING))+
  geom_point(colour = '#1b9e77', size=2, alpha=0.5)+
  theme_classic()+
  labs(y='Kinship Coefficient (Rxy)', x='Distance (m)')

sisobd = ggplot(cj, aes(x=dist,y=KING))+
  geom_point(colour = '#d95f02', size=2, alpha=0.5)+
  theme_classic()+
  labs(y='Kinship Coefficient (Rxy)', x='Distance (m)')


plot_grid(misobd, sisobd)

