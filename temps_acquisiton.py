def temps_acquisition():
    jour_fin = float(input("Veuillez entrer le jour à la fin du temps d'acquisition:"))
    heure_fin = float(input("Veuillez entrer l'heure à la fin du temps d'acquisition:"))
    minute_fin = float(input("Veuillez entrer les minutes à la fin du temps d'acquisition:"))
    jour_debut = float(input("Veuillez entrer le jour au début du temps d'acquisition:"))
    heure_debut = float(input("Veuillez entrer l'heure au début du temps d'acquisition:"))
    minute_debut = float(input("Veuillez entrer les minutes au début du temps d'acquisition:"))

    temps = (jour_fin - jour_debut)*24 + (heure_fin -heure_debut) + (minute_fin - minute_debut)/60.0
    temps_acquisition = round(temps, 2)
    print(f"Le temps d'acquisition est de {temps_acquisition} heures.")
    return temps_acquisition
    
temps_acquisition()