import * as param from "../variables/const.js"

const render = (arr, label, index, title1 = `Заголовок`) => {
    const title = `<h2 style="margin: 25px 0 0 0; font-size: 24px; padding:10px 5px; ">${title1}</h3>`
    const app = document.querySelector('#app')
    const box = document.createElement('div')
    const table = document.createElement('table')
    table.className = "data"
    const header = document.createElement('tr')
    const td = (element) => `<td class="data__element">${element}</td>`
    header.insertAdjacentHTML('beforeend', td(''))
    for (let i = index; i < arr[0].length + index; i++) {
        header.insertAdjacentHTML('beforeend', td(i))
    }
    for (let i = 0; i < arr.length; i++) {
        const tr = document.createElement('tr');
        tr.insertAdjacentHTML('beforeend', td(label[i]))
        arr[i].forEach(element => tr.insertAdjacentHTML('beforeend', td(math.format(element, { notation: 'exponential', precision: 5 }))))
        table.insertAdjacentElement('beforeend', tr)
    }
    table.insertAdjacentElement('afterbegin', header)
    box.insertAdjacentHTML('beforeend', title)
    box.insertAdjacentElement('beforeend', table)
    app.insertAdjacentElement('beforeend', box)
}

const myId = (app) => {
    app.insertAdjacentHTML('afterbegin', `<section class="id">
<h2 class="id__title">Мои ИД:</h2>
<dl class="id__list">
    <div class="id__box">
        <dt class="id__item">Группа:</dt>
        <dd class="id__item">СМ3-72</dd>
    </div>
    <div class="id__box">
        <dt class="id__item">Вариант N:</dt>
        <dd class="id__item">${param.N}</dd>
    </div>
    <div class="id__box">
        <dt class="id__item">h, [м]:</dt>
        <dd class="id__item">${param.h}</dd>
    </div>
    <div class="id__box">
        <dt class="id__item">M, [-]:</dt>
        <dd class="id__item">${param.M_inf}</dd>
    </div>
    <div class="id__box">
        <dt class="id__item">α, [град]:</dt>
        <dd class="id__item">${param.alfa}</dd>
    </div>
    <div class="id__box">
        <dt class="id__item">T, [K]:</dt>
        <dd class="id__item">${param.T_st}</dd>
    </div>
    <div class="id__box">
        <dt class="id__item">Reкр, [-]:</dt>
        <dd class="id__item">${param.Re_k}</dd>
    </div>
    <div class="id__box">
        <dt class="id__item">с, [м]:</dt>
        <dd class="id__item">${param.c}</dd>
    </div>
    <div class="id__box">
        <dt class="id__item">b, [м]:</dt>
        <dd class="id__item">${param.b}</dd>
    </div>
</dl>
</section>`)
}

export { render, myId }